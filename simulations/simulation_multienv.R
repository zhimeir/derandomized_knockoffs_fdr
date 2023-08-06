#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
seedA <- as.integer(args[1])
amp <- as.integer(args[2])

suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(knockoff))
suppressPackageStartupMessages(library(tidyverse))
source("utils.R")
source("MEKF/functions_multienv.R")
source("MEKF/accumulation_test_functions.R")
source("MEKF/importance_stats.R")

## The directory to save the results
save_dir <- sprintf("../results/simulation_multienv")
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
nenv <- 2

## The Y|X model for two environments
nonzero <- list()
beta_true <- list()
y.sample <- list()
nonzero[[1]] <- 1:k
nonzero[[2]] <- 1:(k+10)
for(i in 1:nenv){
  beta_true[[i]] <- amp / 10 * (1:p %in% nonzero[[i]]) / sqrt(n)
  y.sample[[i]] <- function(X) X %*% beta_true[[i]] + rnorm(n)
}

diags <- knockoff::create.solve_asdp(Sigma)

## Generating data
all_res <- data.frame()
set.seed(seedA)
  
Xlist <- list()
Ylist <- list()
for (i in 1:nenv){
  Xlist[[i]] <- matrix(rnorm(n * p),n) %*% chol(Sigma)
  Ylist[[i]] <- y.sample[[i]](Xlist[[i]])
} 
  
## Vanilla  knockoff
Xklist <- list()
for(i in 1:nenv){
  Xklist[[i]] <- create.gaussian(Xlist[[i]], rep(0,p), Sigma, diag_s = diags)
}
  
W_matrix <- compute_stats_with_prior(Xlist, Xklist, Ylist)
W_sign <- 2*((apply(W_matrix, 1, min) > 0) - 0.5)
W_prod <- abs(apply(W_matrix, 1, prod))
W <- W_sign * W_prod
tau <- knockoff.threshold(W, fdr = alpha, offset = 1)
rej <- which(W >= tau)
fdp <- sum(beta_true[[1]][rej]==0) / max(length(rej), 1)
power <- sum(beta_true[[1]][rej]!=0) / k
vkn_res <- data.frame(method = "vanilla", power = power, fdp = fdp, seed = seedA)

## Derandomized Knockoffs
res <- multienv_ekn(Xlist, Ylist, nenv, M, 
           alpha, alpha / 2, rep(0,p), Sigma, 
           diags, family = "gaussian", offset = 1)
rej <- res$rej
fdp <- sum(beta_true[[1]][rej]==0) / max(length(rej), 1)
power <- sum(beta_true[[1]][rej]!=0) / max(k,1)
mkn_res <- data.frame(method = "multiple", power = power, fdp = fdp, seed = seedA)

## Save the results
all_res <- rbind(all_res, vkn_res, mkn_res)


out_dir <- sprintf("%s/res_amp_%d_seedA_%d.csv", save_dir, amp, seedA)
write_csv(all_res, out_dir)

