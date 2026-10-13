#!/usr/bin/env Rscript
# =============================================================================
# 04b - Regional beta-binomial model
#
# Fits a beta-binomial model with each population as an independent observation: regional intervals and tests.
#
# Inputs:  results/freqs_with_ci.tsv, results/biomedical_samples.tsv
# Outputs: results/regional_betabinomial.tsv
# Usage:   Rscript 04b_regional_betabinomial.R [config.yaml]
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(yaml) })

args <- commandArgs(trailingOnly = TRUE)
cfg  <- yaml::read_yaml(if (length(args) >= 1) args[1] else "config.yaml")
P <- cfg$params
region_levels <- unlist(P$region_levels)
level <- P$ci_level

# ---- Model functions ----
LO <- -12; HI <- 12
RLO <- -12; RHI <- 4

bb_ll <- function(c, n, mu, rho) {
  a <- mu * (1 - rho) / rho
  b <- (1 - mu) * (1 - rho) / rho
  sum(lchoose(n, c) + lbeta(c + a, n - c + b) - lbeta(a, b))
}

fit_bb <- function(c, n, g, K, fix = NULL, start = NULL) {
  free <- setdiff(seq_len(K), if (is.null(fix)) integer(0) else fix$k)
  emp  <- vapply(seq_len(K), function(k)
            qlogis((sum(c[g == k]) + 0.5) / (sum(n[g == k]) + 1)), numeric(1))
  nll <- function(par) {
    eta <- numeric(K)
    eta[free] <- par[seq_along(free)]
    if (!is.null(fix)) eta[fix$k] <- fix$eta
    -bb_ll(c, n, plogis(eta[g]), plogis(par[length(par)]))
  }
  lower <- c(rep(LO, length(free)), RLO); upper <- c(rep(HI, length(free)), RHI)
  starts <- if (!is.null(start)) list(start) else
    lapply(c(-6, -3, -1), function(r0) c(pmin(pmax(emp[free], LO), HI), r0))
  best <- NULL
  for (st in starts) {
    o <- optim(st, nll, method = "L-BFGS-B", lower = lower, upper = upper)
    if (is.null(best) || o$value < best$value) best <- o
  }
  eta <- numeric(K)
  eta[free] <- best$par[seq_along(free)]
  if (!is.null(fix)) eta[fix$k] <- fix$eta
  list(ll = -best$value, eta = eta, mu = plogis(eta),
       rho = plogis(best$par[length(best$par)]), lrho = best$par[length(best$par)])
}

profile_ci <- function(c, n, g, K, k, full) {
  crit <- qchisq(level, 1) / 2
  warm <- c(full$eta[-k], full$lrho)
  f <- function(e) full$ll - fit_bb(c, n, g, K, fix = list(k = k, eta = e), start = warm)$ll - crit
  e_hat <- min(max(full$eta[k], LO), HI)
  lo <- if (e_hat <= LO + 1e-6 || f(LO) < 0) 0 else plogis(uniroot(f, c(LO, e_hat), tol = 1e-4)$root)
  hi <- if (e_hat >= HI - 1e-6 || f(HI) < 0) 1 else plogis(uniroot(f, c(e_hat, HI), tol = 1e-4)$root)
  c(lo, hi)
}

lrt_region <- function(c, n, in_r) {
  two <- fit_bb(c, n, ifelse(in_r, 1L, 2L), 2)
  one <- fit_bb(c, n, rep(1L, length(c)), 1)
  pchisq(max(0, 2 * (two$ll - one$ll)), df = 1, lower.tail = FALSE)
}

# ---- Simulation check ----
set.seed(2)
sim_g <- rep(1:3, each = 300); sim_n <- sample(16:64, 900, TRUE)
sim_mu <- c(0.10, 0.30, 0.50)[sim_g]; sim_rho <- 0.05
sim_p <- rbeta(900, sim_mu * (1 - sim_rho) / sim_rho, (1 - sim_mu) * (1 - sim_rho) / sim_rho)
sf <- fit_bb(rbinom(900, sim_n, sim_p), sim_n, sim_g, 3)
cat(sprintf("Sanity check (truth mu 0.10/0.30/0.50, rho 0.050): fitted mu %s, rho %.3f\n",
            paste(sprintf("%.3f", sf$mu), collapse = "/"), sf$rho))

# ---- Data ----
d   <- fread("results/freqs_with_ci.tsv", colClasses = list(character = "group_code"))
map <- unique(fread("results/biomedical_samples.tsv")[, .(CLST_name, Group)])
pop <- merge(d[scale == "clst", .(SNP, pop = group_name, c, n)],
             map, by.x = "pop", by.y = "CLST_name")
pop[, g := match(Group, region_levels)]
stopifnot(!anyNA(pop$g))

# ---- Regional fits ----
res <- rbindlist(lapply(split(pop, pop$SNP), function(x) {
  x <- x[order(g)]
  K <- length(region_levels)
  full <- fit_bb(x$c, x$n, x$g, K)
  polymorphic <- sum(x$c) > 0 && sum(x$c) < sum(x$n)
  rbindlist(lapply(seq_len(K), function(k) {
    ci <- profile_ci(x$c, x$n, x$g, K, k, full)
    in_r <- x$g == k
    data.table(SNP = x$SNP[1], region = region_levels[k], n_pops = sum(in_r),
               rho = full$rho, bb_mean = full$mu[k], bb_lo = ci[1], bb_hi = ci[2],
               p_bb = if (polymorphic) lrt_region(x$c, x$n, in_r) else NA_real_)
  }))
}))
fwrite(res, "results/regional_betabinomial.tsv", sep = "\t", na = "NA")
cat(sprintf("Beta-binomial regional model: %d variants x %d regions; median rho %.3f\n",
            uniqueN(res$SNP), length(region_levels), median(unique(res[, .(SNP, rho)])$rho)))
cat("Wrote results/regional_betabinomial.tsv\n")
