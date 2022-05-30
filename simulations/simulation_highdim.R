#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
seedA<- as.integer(args[1])
amp <- as.integer(args[2])
if(is.na(seedA)) seedA <- 1
if(is.na(amp)) amp <- 10

suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(knockoff))
suppressPackageStartupMessages(library(tidyverse))
source("../utils/utils.R")

## The directory to save the results
save_dir <- sprintf("../results/simulation_highdim")
if(!dir.exists(save_dir)){
  dir.create(save_dir)
}

## Parameters
set.seed(24601)
n <- 300
p <- 500
k <- 50
alpha <- 0.1
rho <- 0.5
M <- 50
Sigma <- toeplitz(rho^(0:(p-1)))
## nonzero <- sample(1:p, k, replace = FALSE)
## beta_true <- amp * (1:p %in% nonzero) * sign(rnorm(p)) / sqrt(n)
nonzero <- seq(1, 2 * k, by = 2)
beta_true <- amp * (1:p %in% nonzero) / sqrt(n)
beta_true[seq(3, 2 * k, 4)] <- -beta_true[seq(3, 2 * k, 4)]

y.sample <- function(X) X %*% beta_true + rnorm(n)
set <- data.frame(vkn = rep(0,p), mkn = rep(0,p), truth = beta_true)
diags <- knockoff::create.solve_asdp(Sigma)
nrep <- 20
all_res <- data.frame()

## Generating data
set.seed(seedA)
X <- matrix(rnorm(n * p),n) %*% chol(Sigma)
Y <- y.sample(X)

## Run the procedures for multiple runs
for(seedB in 1:nrep){
  cat(sprintf("Running the %d-th rep.\n", seedB)) 
  xk_seed <- 100 + (seedA - 1) * nrep + seedB
  set.seed(xk_seed)

  ## Vanilla  knockoff
  Xk <- create.gaussian(X, rep(0,p), Sigma, diag_s = diags)
  W <- stat.glmnet_coefdiff(X, Xk, Y)
  tau <- knockoff.threshold(W, fdr = alpha, offset = 1)
  rej <- which(W >= tau)
  fdp <- sum(beta_true[rej]==0) / max(length(rej), 1)
  power <- sum(beta_true[rej]!=0) / k
  vkn_res <- data.frame(method = "vanilla", power = power, fdp = fdp, seedB = seedB)
  set$vkn[rej] <- set$vkn[rej] + 1

  ## FDR-derandomized knockoffs
  res <- ekn(X, Y, M, alpha, alpha / 2, Sigma, diags, "gaussian", offset = 1)
  rej <- res$rej
  fdp <- sum(beta_true[rej]==0) / max(length(rej), 1)
  power <- sum(beta_true[rej]!=0) / k
  mkn_res <- data.frame(method = "multiple", power = power, fdp = fdp, seedB = seedB)
  set$mkn[rej] <- set$mkn[rej] + 1

  ## Save the results
  all_res <- rbind(all_res, vkn_res, mkn_res)
}

## Save the outcome
out_dir <- sprintf("%s/res_amp_%d_seedA_%d.csv", save_dir, amp, seedA)
write_csv(all_res, out_dir)

out_dir <- sprintf("%s/res_amp_%d_seedA_%d_set.csv", save_dir, amp, seedA)
write_csv(set, out_dir)

