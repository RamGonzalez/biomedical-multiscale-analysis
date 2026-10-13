#!/usr/bin/env python
# =============================================================================
# 02 - Allele frequencies
#
# Counts copies of the effect allele (c) and called chromosomes (n) per population and region with snputils.
#
# Inputs:  results/effect_variant_table.tsv, results/biomedical_samples.tsv, genotypes .bed/.bim/.fam (config.yaml)
# Outputs: results/freq_clst.tsv, results/freq_group.tsv, results/orientation_report.tsv
# Usage:   python 02_compute_freqs_snputils.py [config.yaml]
# =============================================================================

import sys, os
import numpy as np
import pandas as pd
import yaml
import snputils as su

cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config.yaml"
with open(cfg_path) as fh:
    cfg = yaml.safe_load(fh)
IN = cfg["inputs"]
os.makedirs("results", exist_ok=True)

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
def revcomp_base(b):
    return COMP.get(str(b).upper(), None)

# ---- Effect variants ----
eff = pd.read_csv("results/effect_variant_table.tsv", sep="\t", dtype=str)
eff = eff[eff["in_ogvp"] == "TRUE"].copy()
eff["REF"] = eff["REF"].str.upper()
eff["other_allele"] = eff["other_allele"].str.upper()
eff["effect_allele"] = eff["effect_allele"].str.upper()
print(f"[02] effect variants (in OGVP): {len(eff)}", flush=True)

# ---- Genotypes ----
def ensure_tab_fam(src_prefix, workdir="results/work"):
    os.makedirs(workdir, exist_ok=True)
    dst_prefix = os.path.join(workdir, "geno")
    for ext in (".bed", ".bim"):
        link = dst_prefix + ext
        if os.path.islink(link) or os.path.exists(link):
            os.remove(link)
        os.symlink(os.path.abspath(src_prefix + ext), link)
    with open(src_prefix + ".fam") as fh:
        rows = [ln.split() for ln in fh if ln.strip()]
    with open(dst_prefix + ".fam", "w") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return dst_prefix

geno_prefix = ensure_tab_fam(IN["geno_prefix"])
print("[02] reading genotypes ...", flush=True)
snp = su.read_snp(geno_prefix + ".bed")
vid = np.asarray(snp.variants_id).astype(str)
print(f"[02] genotypes: {snp.n_snps} SNPs x {snp.n_samples} samples", flush=True)

# ---- Subset to effect variants ----
pos_of = {}
for i, v in enumerate(vid):
    pos_of.setdefault(v, i)
idx, kept_snp = [], []
for s in eff["SNP"]:
    if s in pos_of:
        idx.append(pos_of[s]); kept_snp.append(s)
    else:
        print(f"[02] WARNING: effect SNP not found in .bim: {s}", flush=True)
idx = np.asarray(idx, dtype=int)
sub = snp.filter_variants(indexes=idx, inplace=False)
sub_id  = np.asarray(sub.variants_id).astype(str)
sub_ref = np.asarray(sub.variants_ref).astype(str)
sub_alt = np.asarray(sub.variants_alt).astype(str)
print(f"[02] subset to {sub.n_snps} effect SNPs", flush=True)

# ---- Sample groups ----
samples = np.asarray(snp.samples).astype(str)
numeric_id = np.array([s.rsplit("_", 1)[-1] for s in samples])

manifest = pd.read_csv("results/biomedical_samples.tsv", sep="\t", dtype=str)
id2clst  = dict(zip(manifest["sample_id"], manifest["CLST"]))
id2group = dict(zip(manifest["sample_id"], manifest["Group"].fillna("UNMAPPED")))
fid   = np.array([id2clst.get(x, "UNMAPPED") for x in numeric_id])
group = np.array([id2group.get(x, "UNMAPPED") for x in numeric_id])
assert not (fid == "UNMAPPED").any(), "samples missing from results/biomedical_samples.tsv"
n_unmapped = int((group == "UNMAPPED").sum())
if n_unmapped:
    print(f"[02] WARNING: {n_unmapped} samples had no Group mapping "
          f"(labelled UNMAPPED, excluded from panel b)", flush=True)

# ---- Allele counts per group ----
def counts_by(labels):
    alt_df, called_df = sub.allele_counts(sample_labels=labels,
                                          return_called=True, as_dataframe=True)
    alt_df = alt_df.reset_index(drop=True)
    called_df = called_df.reset_index(drop=True)
    return alt_df, called_df

alt_a, called_a = counts_by(fid)
alt_b, called_b = counts_by(group)

# ---- Allele orientation ----
eff_by_snp = eff.set_index("SNP")
orient = []
for i, s in enumerate(sub_id):
    e_ref = eff_by_snp.at[s, "other_allele"]
    e_alt = eff_by_snp.at[s, "effect_allele"]
    s_ref, s_alt = sub_ref[i].upper(), sub_alt[i].upper()
    if s_alt == e_alt and s_ref == e_ref:
        status, flip = "aligned", False
    elif s_alt == e_ref and s_ref == e_alt:
        status, flip = "flipped", True
    elif revcomp_base(s_alt) == e_alt and revcomp_base(s_ref) == e_ref:
        status, flip = "strand_aligned", False
    elif revcomp_base(s_alt) == e_ref and revcomp_base(s_ref) == e_alt:
        status, flip = "strand_flipped", True
    else:
        status, flip = "MISMATCH", None
    orient.append((s, e_ref, e_alt, s_ref, s_alt, status, flip))

orient_df = pd.DataFrame(orient, columns=[
    "SNP", "effect_ref", "effect_allele", "snputils_ref", "snputils_alt",
    "orientation_status", "flip"])
orient_df.to_csv("results/orientation_report.tsv", sep="\t", index=False)
bad = orient_df[orient_df["orientation_status"] == "MISMATCH"]
print(f"[02] orientation: "
      + " | ".join(f"{k}={v}" for k, v in orient_df['orientation_status'].value_counts().items()),
      flush=True)
if len(bad):
    print(f"[02] *** {len(bad)} variant(s) with UNRESOLVABLE alleles — quarantined:", flush=True)
    print(bad.to_string(index=False), flush=True)

flip_map = dict(zip(orient_df["SNP"], orient_df["flip"]))
status_map = dict(zip(orient_df["SNP"], orient_df["orientation_status"]))

# ---- Frequency tables ----
meta_cols = ["SNP", "rsid", "REF", "effect_allele", "category", "source", "snp_label"]
meta = eff[meta_cols].copy()

def build_long(alt_df, called_df, scale_name):
    alt_df = alt_df.copy();   alt_df.insert(0, "SNP", sub_id)
    called_df = called_df.copy(); called_df.insert(0, "SNP", sub_id)
    a_long = alt_df.melt(id_vars="SNP", var_name="group_code", value_name="alt_count")
    n_long = called_df.melt(id_vars="SNP", var_name="group_code", value_name="n")
    df = a_long.merge(n_long, on=["SNP", "group_code"])
    df["group_code"] = df["group_code"].astype(str)
    df = df[df["group_code"] != "UNMAPPED"]
    df["flip"] = df["SNP"].map(flip_map)
    df["orientation_status"] = df["SNP"].map(status_map)
    df = df[df["orientation_status"] != "MISMATCH"]
    df["c"] = np.where(df["flip"], df["n"] - df["alt_count"], df["alt_count"]).astype("float")
    df["n"] = df["n"].astype("float")
    df["effect_freq"] = np.where(df["n"] > 0, df["c"] / df["n"], np.nan)
    df["scale"] = scale_name
    df = df.merge(meta, on="SNP", how="left")
    return df[["scale", "SNP", "rsid", "REF", "effect_allele", "category", "source",
               "snp_label", "group_code", "c", "n", "effect_freq", "orientation_status"]]

long_a = build_long(alt_a, called_a, "clst")
long_b = build_long(alt_b, called_b, "group")

clst_map = dict(zip(manifest["CLST"], manifest["CLST_name"]))
long_a["group_name"] = long_a["group_code"].map(clst_map).fillna(long_a["group_code"])
long_b["group_name"] = long_b["group_code"]

long_a.to_csv("results/freq_clst.tsv", sep="\t", index=False, na_rep="NA")
long_b.to_csv("results/freq_group.tsv", sep="\t", index=False, na_rep="NA")
print(f"[02] wrote results/freq_clst.tsv  ({long_a['SNP'].nunique()} SNPs x "
      f"{long_a['group_code'].nunique()} CLST)", flush=True)
print(f"[02] wrote results/freq_group.tsv ({long_b['SNP'].nunique()} SNPs x "
      f"{long_b['group_code'].nunique()} groups)", flush=True)
