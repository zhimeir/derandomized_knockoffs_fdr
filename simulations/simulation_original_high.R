#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
seedA<- as.integer(args[1])
amp <- as.integer(args[2])
if(is.na(seedA)) seedA <- 1
if(is.na(amp)) amp <- 30

suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(knockoff))
suppressPackageStartupMessages(library(tidyverse))
source("../utils/utils.R")

## The directory to save the results
save_dir <- sprintf("../results/simulation_original_high")
if(!dir.exists(save_dir)){
  dir.create(save_dir)
}

## Parameters
set.seed(24601)
n <- 3000
p <- 6000
k <- 60
alpha <- 0.1
rho <- 0
M <- 50
mu <- rep(0,p)
Sigma <- toeplitz(rho^(0:(p-1)))
nonzero <- sample(1:p, k, replace = FALSE)
beta_true <- amp * sign(rnorm(p)) * (1:p %in% nonzero) / sqrt(n)

y.sample <- function(X) X %*% beta_true + rnorm(n)
diags <- knockoff::create.solve_asdp(Sigma)
all_res <- data.frame()

## Generating data
set.seed(seedA)
X <- matrix(rnorm(n * p),n) %*% chol(Sigma)
Y <- y.sample(X)


## Vanilla  knockoff
Xk <- create.gaussian(X, mu, Sigma, diag_s = diags)
W <- stat.glmnet_coefdiff(X, Xk, Y)
tau <- knockoff.threshold(W, fdr = alpha, offset = 1)
rej <- which(W >= tau)
fdp <- sum(beta_true[rej]==0) / max(length(rej), 1)
power <- sum(beta_true[rej]!=0)/max(k,1)
all_res <- rbind(all_res,
  data.frame(method = "vanilla", power = power, fdp = fdp))

## Compute knockoff e-values 
res <- ekn(X, Y, M, alpha/2, mu, Sigma, diags, "gaussian", offset = 1)
E <- res$E
rej <- ebh(E, alpha)$rej
fdp <- sum(beta_true[rej]==0) / max(length(rej), 1)
power <- sum(beta_true[rej]!=0) / k
all_res <- rbind(all_res,
    data.frame(method = "multiple", power = power, fdp = fdp))

out_dir <- sprintf("%s/res_amp_%d_seedA_%d.csv", save_dir, amp, seedA)
write_csv(all_res, out_dir)

