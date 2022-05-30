#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
seedA <- as.integer(args[1])
amp <- as.numeric(args[2])

suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(gam))
suppressPackageStartupMessages(library(knockoff))
suppressPackageStartupMessages(library(adaptiveKnockoff))
suppressPackageStartupMessages(library(tidyverse))
source("../utils/utils.R")

## The directory to save the results
save_dir <- sprintf("../results/simulation_sideinfo")
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
Sigma <- toeplitz(rho^(0:(p-1)))
nonzero <- 1:k
beta_true <- amp * (1:p %in% nonzero) * sign(rnorm(p)) / 2 / sqrt(n)
y.sample <- function(X) X %*% beta_true + rnorm(n)
diags <- knockoff::create.solve_asdp(Sigma)
set <- data.frame(vkn = rep(0,p), wmkn = rep(0,p), mkn = rep(0,p), truth = beta_true)
nrep <- 1

## Generating data
all_res <- data.frame()

for(seedB in 1:nrep){
  
  seed <- (seedA - 1) * nrep + seedB
  set.seed(seed)

  X <- matrix(rnorm(n * p),n) %*% chol(Sigma)
  Y <- y.sample(X)
  
  ## Vanilla  knockoff
  Xk <- create.gaussian(X, rep(0,p), Sigma, diag_s = diags)
  W <- stat.glmnet_coefdiff(X, Xk, Y, family = "gaussian")
  tau <- knockoff.threshold(W, fdr = alpha, offset = 1)
  rej <- which(W >= tau)
  fdp <- sum(beta_true[rej]==0) / max(length(rej), 1)
  power <- sum(beta_true[rej]!=0) / k
  vkn_res <- data.frame(method = "vanilla", power = power, fdp = fdp, seed = seed)
  set$vkn[rej] <- set$vkn[rej] + 1


  ## Weighted derandomized knockoffs
  weight <- exp(-(1:p))
  res <- weighted_ekn(X, Y, M, alpha, alpha/2, weight, Sigma, diags)
  rej <- res$rej
  fdp <- sum(beta_true[rej]==0) / max(length(rej), 1)
  power <- sum(beta_true[rej]!=0) / k
  wmkn_res <- data.frame(method = "weighted_multiple", 
                         power = power, fdp = fdp, seed = seed)
  set$wmkn[rej] <- set$wmkn[rej] + 1


  ## Derandomized knockoffs
  res <- ekn(X, Y, M, alpha, alpha / 2, Sigma, diags, 
             family = "gaussian", offset = 1)
  rej <- res$rej
  fdp <- sum(beta_true[rej]==0) / max(length(rej), 1)
  power <- sum(beta_true[rej]!=0) / k
  mkn_res <- data.frame(method = "multiple", 
                         power = power, fdp = fdp, seed = seed)
  set$mkn[rej] <- set$mkn[rej] + 1

  ## Save the results
  all_res <- rbind(all_res, vkn_res, mkn_res, wmkn_res)
}


out_dir <- sprintf("%s/res_amp_%.1f_seedA_%d.csv", save_dir, amp, seedA)
write_csv(all_res, out_dir)

out_dir <- sprintf("%s/res_amp_%.1f_seedA_%d_set.csv", save_dir, amp, seedA)
write_csv(set, out_dir)
