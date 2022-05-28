ekn <- function(X, Y, M, alpha, gamma, Sigma, diags, 
                family = "gaussian", offset = 1){

  n <- dim(X)[1]
  p <- dim(X)[2]

  E <- matrix(0, M, p)

  for(m in 1:M){
    Xk <- create.gaussian(X, rep(0,p), Sigma, diag_s = diags)
    W <- stat.glmnet_coefdiff(X, Xk, Y, family = family)

    ## Test different threshold
    tau <- knockoff.threshold(W, fdr =  gamma, offset = offset)
    E[m,] <- (W >= tau) / (1 + sum(W <= -tau))

  }
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

  return(list(rej = rej, E = E))

}


boost_ekn <- function(X, Y, M, alpha, gamma, diags, offset = 0, one_var = TRUE){

  n <- dim(X)[1]
  p <- dim(X)[2]

  E <- matrix(0, M, p)

  for(m in 1:M){

    Xk <- create.gaussian(X, rep(0,p), Sigma, diag_s = diags)
    mdl <- cv.glmnet(cbind(X, Xk), Y, alpha=1)
    cvlambda <- mdl$lambda.min
    beta_hat <- mdl$glmnet.fit$beta[,mdl$lambda ==mdl$lambda.min]
    T <- abs(beta_hat[1:p])
    T_tilde <- abs(beta_hat[(p+1):(2*p)])
    W <- T-T_tilde

    ## Test different threshold
    tau <- knockoff.threshold(W, fdr =  gamma, offset = offset)
    E[m,] <- (W >= tau) / (1 + sum(W <= -tau))

  }

  E <- colMeans(E)
  if(one_var){
    U <- runif(1)
  }else{
    U <- runif(p)
  }
  E <- E / U
  E_ord <- order(E, decreasing = TRUE)
  E <- sort(E, decreasing = TRUE)
  comp <- E >= (1 / alpha / (1:p))
  id <- max(which(comp > 0))
  if(id > 0){
    rej <- E_ord[1:id]
  }else{
    rej <- NULL
  }

  return(list(rej = rej, E = E))

}

compute_mmi_s <- function(Sigma){
  p <- dim(Sigma)[1]
  ones <- rep(1,p)
  s <- Variable(p)
  
    obj <- Maximize(t(ones)%*%log(s)+log_det(2*Sigma-diag(s)))
    constraints <- list(
      s>=0,
      2*Sigma-diag(s)>=0
    )
    prob <- Problem(obj, constraints)
    res <- psolve(prob,solver = "SCS")
    s <- as.vector(res$getValue(s))
  return(s)  
}


weighted_ekn <- function(X, Y, M, alpha, gamma, weight, 
                         Sigma, diags, offset = 1, family = "gaussian"){

  n <- dim(X)[1]
  p <- dim(X)[2]

  E <- matrix(0, M, p)

  for(m in 1:M){
    Xk <- create.gaussian(X, rep(0,p), Sigma, diag_s = diags)
    W <- stat.glmnet_coefdiff(X, Xk, Y, family = family)        

    ## Test different threshold
    tau <- knockoff.threshold(W, fdr =  gamma, offset = offset)
    E[m,] <- (W >= tau) * weight / (weight + sum((W <= -tau) * weight))
  }

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

  return(list(rej = rej, E = E))
}


adaptive_ekn <- function(X, Y, M, alpha, gamma, U, 
                         Sigma, diags, family = "gaussian",
                         offset = 1){

  n <- dim(X)[1]
  p <- dim(X)[2]

  E <- matrix(0, M, p)

  for(m in 1:M){
    Xk <- create.gaussian(X, rep(0,p), Sigma, diag_s = diags)
    W <- stat.glmnet_coefdiff(X, Xk, Y, family = family)    

    ## Test different threshold
    akn_res <- filter_gam(W, U, gamma)
    if_remain <- (1:p) %in% akn_res$unrevealed_id
    E[m,] <- if_remain * (W > 0) / (1 + sum((W < 0) * if_remain))
  }

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

  return(list(rej = rej, E = E))

}

pfer_kn <- function(X, Y, v0, M, tau, Sigma, diags, family = "gaussian"){

  n <- nrow(X)
  p <- ncol(X)
  
  freq <- rep(0, p)
  for (m in 1:M){
    v <- floor(v0)+rbinom(1,1,v0-floor(v0))
    Xk <- create.gaussian(X,rep(0,p),Sigma,diag_s = diags)
    W <- stat.glmnet_coefdiff(X, Xk, Y, family = family)
    order_w <- order(abs(W),decreasing = TRUE)
    sorted_w <- W[order_w]
    negid <- which(sorted_w<0)
    if(v>0){
      if(sum(W<0)<v){
        S <- which(W>0)
        freq[S] <- freq[S]+1
      }else{
        TT <- negid[v]
        S <- which(sorted_w[1:TT]>0)
        S <- order_w[S]
        freq[S] <- freq[S]+1
        }
      }
  }
  freq <- freq/M
  S <- which(freq>=tau)

return(list(rej = S, freq = freq))
}


multienv_ekn <- function(Xlist, Ylist, nenv, 
                         M, alpha, gamma, Sigma, diags, 
                         family = "gaussian", offset = 1){

  n <- dim(Xlist[[1]])[1]
  p <- dim(Xlist[[1]])[2]

  E <- matrix(0, M, p)

  for(m in 1:M){
    
    Xklist <- list()
    for(i in 1:nenv){
      Xklist[[i]] <- create.gaussian(Xlist[[i]], rep(0,p), Sigma, diag_s = diags)
    }
    W_matrix <- compute_stats_with_prior(Xlist, Xklist, Ylist)
    W_sign <- 2*((apply(W_matrix, 1, min) > 0) - 0.5)
    W_prod <- abs(apply(W_matrix, 1, prod))
    W <- W_sign * W_prod

    ## Test different threshold
    tau <- knockoff.threshold(W, fdr =  gamma, offset = offset)
    E[m,] <- (W >= tau) / (1 + sum(W <= -tau))

  }
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

  return(list(rej = rej, E = E))

}
