#!/usr/bin/env Rscript
#
# Author: Zhimei Ren
# Apply derandomized knockoffs to the HIV resistance dataset
#
# Latest update: 2022-05-15
#

#########################
### process arguments ###
#########################
args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1])
set.seed(seed) 

## other parameters
class <- "PI"
drug <- "LPV"
alpha <- 0.1
M <- 100         # number of knockoffs

######################
### load libraries ###
######################
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(knockoff))
suppressPackageStartupMessages(library(glmnet))

#################
### load data ###
#################
data_file <- sprintf("../data/%s_%s.txt", class, drug)
data <- read_delim(data_file, delim = " ", col_types = cols())
Y <- data$Y
X <- data[,-1]
p <- ncol(X)


###################
### normalize X ###
###################
standardization <- function(x){
  x <- x - mean(x)
  x <- x / sd(x)
  return(x)
}
X <- apply(X, 2, standardization)

############################
### generating knockoffs ###
############################
E <- matrix(0, M, p)
all_rej <- c()
for(m in 1:M){

  ## generate knockoffs and compute feature importance statistics
  Xk <- create.second_order(as.matrix(X), shrink = FALSE)
  W <- stat.glmnet_coefdiff(X = X, X_k = Xk, y = Y)

  ## derandomized knockoffs
  tau <- knockoff.threshold(W, fdr = alpha / 2, offset = 1)
  E[m,] <- (W >= tau) / (1 + sum(W <= -tau))

  ## vanilla knockoffs
  tau <- knockoff.threshold(W, fdr = alpha, offset = 1)
  
  discoveries <- data.frame(rej = which(W >= tau))
  out_dir <- sprintf("../results/vanilla_%d.txt", m)
  write_delim(discoveries, out_dir, delim = " ")

  if(m %% 20 == 1){
    cat(sprintf("Implementing the %d-th run.\n", m))
  }
}

## run e-BH
E <- colMeans(E)
E_ord <- order(E, decreasing = TRUE)
E <- sort(E, decreasing = TRUE)
comp <- E >= (1 / alpha / (1:p))
id <- max(which(comp > 0))
if(id > 0){
  rej <- E_ord[1:id]
}else{
  rej <- NULL
}


rej <- data.frame(rej = rej)
out_dir <- sprintf("../results/multiple_%d.txt", seed)
write_delim(rej, out_dir, delim = " ")


