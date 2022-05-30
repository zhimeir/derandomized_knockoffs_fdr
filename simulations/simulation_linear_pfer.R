#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
seedA<- as.integer(args[1])
amp <- as.numeric(args[2])
v0 <- as.numeric(args[3])

suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(knockoff))
suppressPackageStartupMessages(library(tidyverse))
source("../utils/utils.R")

## The directory to save the results
save_dir <- sprintf("../results/simulation_linear_pfer")
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
nonzero <- seq(1, p, by = 2)
beta_true <- amp * (1:p %in% nonzero) / sqrt(n)
beta_true[seq(3, p ,by = 4)] = - beta_true[seq(3, p, by = 4)]
y.sample <- function(X) X %*% beta_true + rnorm(n)
diags <- knockoff::create.solve_asdp(Sigma)
all_res <- data.frame()
nrep <- 2

tau <- 0.5

## Run the procedures for multiple runs
for(seedB in 1:nrep){
  
  seed <- (seedA - 1) * nrep + seedB
  set.seed(seed)

  X <- matrix(rnorm(n * p),n) %*% chol(Sigma)
  Y <- y.sample(X)
  
  ## pfer knockoffs
  res <- pfer_kn(X, Y, v0, M, tau, Sigma, diags, "gaussian")
  rej <- res$rej
  fdp <- sum(beta_true[rej]==0) / max(length(rej), 1)
  power <- sum(beta_true[rej]!=0) / k
  vkn_res <- data.frame(method = "pfer", power = power, fdp = fdp, seed = seed)

  ## Save the results
  all_res <- rbind(all_res, vkn_res)
}

out_dir <- sprintf("%s/res_amp_%.1f_seedA_%d_v_%.2f.csv", save_dir, amp, seedA, v0)
write_csv(all_res, out_dir)


