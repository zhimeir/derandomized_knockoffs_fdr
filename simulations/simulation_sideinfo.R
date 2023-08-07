#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
seedA <- as.integer(args[1])
amp <- as.integer(args[2])

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
mu <- rep(0,p)
Sigma <- toeplitz(rho^(0:(p-1)))
nonzero <- 1:k
beta_true <- amp / 10 * (1:p %in% nonzero) / sqrt(n)
y.sample <- function(X) X %*% beta_true + rnorm(n)
diags <- knockoff::create.solve_asdp(Sigma)

## Generating data
set.seed(seedA)
all_res <- data.frame()

X <- matrix(rnorm(n * p),n) %*% chol(Sigma)
Y <- y.sample(X)
  
## Vanilla  knockoff
Xk <- create.gaussian(X, mu, Sigma, diag_s = diags)
W <- stat.glmnet_coefdiff(X, Xk, Y, family = "gaussian")
tau <- knockoff.threshold(W, fdr = alpha, offset = 1)
rej <- which(W >= tau)
fdp <- sum(beta_true[rej]==0) / max(length(rej), 1)
power <- sum(beta_true[rej]!=0) / k
vkn_res <- data.frame(method = "vanilla", power = power, fdp = fdp, seed = seedA)

## Adaptive Knockoffs
res <- filter_gam(W, 1:p, alpha)
rej <- res$rejs[[1]]
fdp <- sum(beta_true[rej]==0) / max(length(rej), 1)
power <- sum(beta_true[rej]!=0) / k
akn_res <- data.frame(method = "adaptive", power = power, fdp = fdp, seed = seedA)

## Derandomized knockoffs
weight <- exp(-(1:p))
res <- weighted_ekn(X, Y, M, alpha, alpha/2, 
                    weight, mu, Sigma, diags)
rej <- res$rej
fdp <- sum(beta_true[rej]==0) / max(length(rej), 1)
power <- sum(beta_true[rej]!=0) / k
wmkn_res <- data.frame(method = "weighted_multiple", 
                         power = power, fdp = fdp, seed = seedA)


## Derandomized knockoffs
E <- ekn(X, Y, M, alpha / 2, mu, Sigma, diags, 
           family = "gaussian", offset = 1)$E
rej <- ebh(E, alpha)$rej
fdp <- sum(beta_true[rej]==0) / max(length(rej), 1)
power <- sum(beta_true[rej]!=0) / k
mkn_res <- data.frame(method = "multiple", 
                       power = power, fdp = fdp, seed = seedA)

## Save the results
all_res <- rbind(all_res, vkn_res, akn_res, mkn_res, wmkn_res)


out_dir <- sprintf("%s/res_amp_%d_seedA_%d.csv", save_dir, amp, seedA)
write_csv(all_res, out_dir)

