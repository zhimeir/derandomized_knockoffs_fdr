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
save_dir <- sprintf("../results/simulation_pfer")
if(!dir.exists(save_dir)){
  dir.create(save_dir)
}

## Parameters
set.seed(24601)
n <- 1000
p <- 800
k <- 80
alpha <- 0.1
rho <- 0.5
M <- 50
mu <- rep(0,p)
Sigma <- toeplitz(rho^(0:(p-1)))
nonzero <- seq(10, p, by = 10)
sign_loc <- seq(20, p, by = 20)
beta_true <- rep(0, p)
beta_true[nonzero] <- rnorm(k, amp / 10, 1) / sqrt(n)
beta_true[sign_loc] <- -beta_true[sign_loc] 

y.sample <- function(X) X %*% beta_true + rnorm(n)
diags <- knockoff::create.solve_asdp(Sigma)
all_res <- data.frame()

## Generating data
set.seed(seedA)
X <- matrix(rnorm(n * p),n) %*% chol(Sigma)
Y <- y.sample(X)

## Run the procedures for multiple runs
stat_mat <- kn_stat(X,Y,M,mu,Sigma,diags)

## FDR-controlling procedure
evals <- kn_evals(stat_mat, alpha / 2, 1)
rej <- ebh(evals, alpha)$rej
fdp <- sum(beta_true[rej]==0) / max(length(rej), 1)
power <- sum(beta_true[rej]!=0) / k
all_res <- rbind(all_res,
                 data.frame(method = "mkn", 
                            power = power, fdp = fdp, v = 0))



## PFER-controlling procedure
v_list <- seq(0,5,by=0.25)
for(v0 in v_list){
  rej <- pfer(stat_mat, v0, 0.5)$rej 
  fdp <- sum(beta_true[rej]==0) / max(length(rej), 1)
  power <- sum(beta_true[rej]!=0) / k
  all_res <- rbind(all_res,
                 data.frame(method = "pkn", power = power, fdp = fdp, v = v0))

}



out_dir <- sprintf("%s/res_amp_%d_seedA_%d.csv", save_dir, amp, seedA)
write_csv(all_res, out_dir)

