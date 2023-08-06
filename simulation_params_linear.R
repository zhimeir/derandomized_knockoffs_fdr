#!/usr/bin/env Rscript

## Process the input 
args <- commandArgs(trailingOnly = TRUE)
seedA <- as.integer(args[1])
amp <- as.integer(args[2])
if(is.na(seedA)) seedA <- 1
if(is.na(amp)) amp <- 30

## Libraries and utility functions
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(knockoff))
suppressPackageStartupMessages(library(tidyverse))
source("utils.R")

## The directory to save the results
save_dir <- sprintf("../results/simulation_alpha")
if(!dir.exists(save_dir)){
  dir.create(save_dir)
}

## Setting up the model
set.seed(24601)
n <- 1000
p <- 800
k <- 80
rho <- 0.5
mu <- rep(0,p)
Sigma <- toeplitz(rho^(0:(p-1)))
nonzero <- seq(10, p, by = 10)
sign_loc <- seq(20, p, by = 20)
beta_true <- rep(0, p)
beta_true[nonzero] <- rnorm(k, amp / 10, 1) / sqrt(n)
beta_true[sign_loc] <- -beta_true[sign_loc] 
y.sample <- function(X) X %*% beta_true + rnorm(n)
diags <- knockoff::create.solve_asdp(Sigma)

## Parameters
M <- 50
alpha <- 0.1
gamma_list <- seq(0.01,alpha*2,by = 0.01)
offset_list <- seq(0,2,by=0.1)

## Generating data
set.seed(seedA)
X <- matrix(rnorm(n * p),n) %*% chol(Sigma)
Y <- y.sample(X)

## Initialization
all_res <- data.frame()

## Vanilla  knockoffs
Xk <- create.gaussian(X, mu, Sigma, diag_s = diags)
W <- stat.glmnet_coefdiff(X, Xk, Y)
tau <- knockoff.threshold(W, fdr = alpha, offset = 1)
rej <- which(W >= tau)
res <- gen_res(rej, beta_true)
all_res <- rbind(all_res, data.frame(method = "vanilla", 
                                     power = res$power, fdp = res$fdp, 
                                     gamma = alpha, offset = 1))

## Compute knockoff statistics 
stat_mat <- kn_stat(X, Y, M, mu, Sigma, diags, "gaussian")

## Loop through the choices of gamma
for(gamma in gamma_list){
  for(offset in offset_list){

    ## rejection set with the naive stopping time
    E <- kn_evals(stat_mat, gamma = gamma, offset = offset, early = FALSE)
    rej <- ebh(E,alpha)$rej
    res <- gen_res(rej, beta_true) 
    all_res <- rbind(all_res, data.frame(method = "multiple", 
                                         power = res$power, fdp = res$fdp, 
                                         gamma = gamma, offset = offset))


    ## rejection set with the early stopping time
    E <- kn_evals(stat_mat, gamma = gamma, offset = offset, early = TRUE)
    rej <- ebh(E,alpha)$rej
    res <- gen_res(rej, beta_true) 
    all_res <- rbind(all_res, data.frame(method = "early", power = res$power, 
                                         fdp = res$fdp, gamma = gamma, offset = offset))
  }  
}


out_dir <- sprintf("%s/res_seedA_%d_amp_%d.csv", save_dir, seedA, amp)
write_csv(all_res, out_dir)


