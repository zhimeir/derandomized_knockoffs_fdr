#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
seedA<- as.integer(args[1])
amp <- as.integer(args[2])
n_ex <- as.integer(args[3])
if(is.na(seedA)) seedA <- 1
if(is.na(amp)) amp <- 5
if(is.na(n_ex)) n_ex <- 500

suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(knockoff))
suppressPackageStartupMessages(library(tidyverse))
source("../utils/utils.R")

## The directory to save the results
save_dir <- sprintf("../results/simulation_robustness")
if(!dir.exists(save_dir)){
  dir.create(save_dir)
}

## Parameters
set.seed(24601)
n <- 600
p <- 100
k <- 50
alpha <- 0.1
rho <- 0.5
M <- 50
mu <- rep(0,p)
Sigma <- toeplitz(rho^(0:(p-1)))
nonzero <- seq(1, p, by = 2)
beta_true <- amp / 10 * (1:p %in% nonzero) / sqrt(n)
beta_true[seq(3, p, by = 4)] <- -beta_true[seq(3, p, by = 4)]

y.sample <- function(X) X %*% beta_true + rnorm(n)
nrep <- 5
all_res <- data.frame()

## Estimate covariance matrix using extra unsupervised data
X_ex <- matrix(rnorm(n_ex * p), n_ex) %*% chol(Sigma)
Sigma_est <- cov(X_ex)
diags <- knockoff::create.solve_asdp(Sigma_est)

## Generating data

## Run the procedures for multiple runs
for(seedB in 1:nrep){
  seed <- (seedA - 1) * nrep + seedB
  set.seed(seed)

  X <- matrix(rnorm(n * p),n) %*% chol(Sigma)
  Y <- y.sample(X)

  ## Vanilla  knockoff
  Xk <- create.gaussian(X, mu, Sigma_est, diag_s = diags)
  W <- stat.glmnet_coefdiff(X, Xk, Y)
  tau <- knockoff.threshold(W, fdr = alpha, offset = 1)
  rej <- which(W >= tau)
  fdp <- sum(beta_true[rej]==0) / max(length(rej), 1)
  power <- sum(beta_true[rej]!=0) / k
  all_res <- rbind(all_res,
                   data.frame(method = "vanilla", power = power,
                              fdp = fdp, seed = seed))

  ## Compute knockoff e-values 
  res <- ekn(X, Y, M, alpha/2, mu, Sigma_est, diags, "gaussian", offset = 1) 
  E <- res$E

  ## rejection set with the naive stopping time
  rej <- ebh(E, alpha)$rej
  fdp <- sum(beta_true[rej]==0) / max(length(rej), 1)
  power <- sum(beta_true[rej]!=0) / k
  all_res <- rbind(all_res,
    data.frame(method = "multiple", power = power, fdp = fdp, seed = seed))
}

out_dir <- sprintf("%s/res_amp_%d_seedA_%d_n_%d.csv", save_dir, amp, seedA, n_ex)
write_csv(all_res, out_dir)

