################################################################
## These functions compute the empirical cross-prior statistics
################################################################
library(glmnet)

## Computes empirical cross-prior statistics 
## Inputs: 
##    Xs: a list of matrix X from different environments
##    X_ks:  a list of knockoffs matrix X_k from different environments
##    ys:  a list of vector y from different environments
## Outputs:
##    A matrix of empirical cross-prior statistics
compute_stats_with_prior = function(Xs, X_ks, ys, family = "gaussian", verbose=FALSE, dfmax=500) {
    n = dim(Xs[[1]])[1]
    p = dim(Xs[[1]])[2]
    num_env = length(Xs)
    W_new_matrix = matrix(rep(0,num_env*p), ncol = num_env)

    for (k in 1:num_env){
        if(verbose) cat(sprintf("Computing statistics for environment %d...\n", k))

        X_stack_mask = NULL
        X_stack_k_mask = NULL
        y_stack = NULL
        env_membership = NULL
        for (k2 in 1:num_env){
            ## Swap original-knockoffs in all environments
            M.swap = matrix(rbinom(n*p,1,1/2), n)
            X_mask = Xs[[k2]] * (1-M.swap) + X_ks[[k2]] * M.swap
            X_k_mask = Xs[[k2]] * M.swap + X_ks[[k2]] * (1-M.swap)
            X_stack_mask = rbind(X_stack_mask, X_mask)
            X_stack_k_mask = rbind(X_stack_k_mask, X_k_mask)
            y_stack = c(y_stack, ys[[k2]])
            env_membership = c(env_membership, rep(k2, dim(Xs[[k2]])[1]))
        }

        ## Fit the lasso on all environments
        X_Xk_stack_mask = cbind(X_stack_mask, X_stack_k_mask)
        ##cv_fit = cv.glmnet(X_Xk_stack_mask[env_membership!=k,], y_stack[env_membership!=k], alpha=0)
        cv_fit = cv.glmnet(X_Xk_stack_mask, y_stack, alpha=0, family = family)
        beta.hat.prior = coef(cv_fit, s="lambda.min")[-1]
        beta.hat.prior = abs(beta.hat.prior[1:p])+abs(beta.hat.prior[(p+1):(2*p)])

        ## Fit the lasso on the data from environment k, tuning weight of prior
        X_env = cbind(Xs[[k]], X_ks[[k]])
        eval_gamma = function(gamma) {
            penalty = 1*(1-gamma)  + gamma * 1 / (0.05+beta.hat.prior)
            cv_fit = cv.glmnet(X_env, ys[[k]], penalty.factor=rep(penalty,2), dfmax=dfmax, family = family)
            idx.min = which.min(cv_fit$cvlo+cv_fit$cvup)
            err = c(cv_fit$cvlo[idx.min], cv_fit$cvup[idx.min])
            return(err)
        }        
        gamma.seq = seq(0,1,length.out=10) # Sequence of prior weights to consider
        err.seq = sapply(gamma.seq, function(gamma) eval_gamma(gamma))

        if(FALSE) {
            tibble(gamma=gamma.seq, low=err.seq[1,], up=err.seq[2,]) %>%
                mutate(mid = (low+up)/2) %>%
                gather(low, up, mid, key="Limit", value="Err") %>%            
                ggplot(aes(x=gamma, y=Err, color=Limit)) +
                geom_point() +
                geom_line()
        }        
        
        ## Find optimal prior weight
        idx.best = which.min(colMeans(err.seq))
        ## Make sure this is significantly better than gamma=0
        if(err.seq[2,idx.best] >= mean(err.seq[,1])) {
            idx.best  = 1
        }
        gamma = gamma.seq[idx.best]
        if(verbose) cat(sprintf("Optimal prior weight: %.3f\n", gamma))
        
        ## Re-fit the lasso on environment k with optimally tuned prior
        penalty = 1*(1-gamma)  + gamma * 1 / (0.05+beta.hat.prior)
        cv_fit = cv.glmnet(X_env, ys[[k]], penalty.factor=rep(penalty,2), dfmax=dfmax, family = family)

        ## Extract importance measures
        beta.hat = coef(cv_fit, s="lambda.min")[-1][1:(2*p)]
        W1 = abs(beta.hat[1:p]) - abs(beta.hat[(p+1):(2*p)])
        W_new_matrix[,k] = W1
    }
    return(W_new_matrix)
}

## Computes empirical cross-prior statistics for transfer learning
## The first environment is treated as the target environment
## Inputs: 
##    Xs: a list of matrix X from different environments
##    X_ks:  a list of knockoffs matrix X_k from different environments
##    ys:  a list of vector y from different environments
## Outputs:
##    A vector of empirical cross-prior statistics
compute_transfer_stats_with_prior = function(Xs, X_ks, ys, family = "gaussian", verbose=FALSE, dfmax=500) {
    ## assume target environment is environment 1
    n = dim(Xs[[1]])[1]
    p = dim(Xs[[1]])[2]
    num_env = length(Xs)

    X_stack = NULL
    X_stack_k = NULL
    y_stack = NULL
    env_membership = NULL
    for (k2 in 2:num_env){
        X_stack = rbind(X_stack, Xs[[k2]])
        X_stack_k = rbind(X_stack_k, X_ks[[k2]])
        y_stack = c(y_stack, ys[[k2]])
        env_membership = c(env_membership, rep(k2, dim(Xs[[k2]])[1]))
    }
    
    ## Fit the lasso on all other environments
    X_Xk_stack = cbind(X_stack, X_stack_k)
    ##cv_fit = cv.glmnet(X_Xk_stack_mask[env_membership!=k,], y_stack[env_membership!=k], alpha=0)
    cv_fit = cv.glmnet(X_Xk_stack, y_stack, alpha=0, family = family)
    beta.hat.prior = coef(cv_fit, s="lambda.min")[-1]
    beta.hat.prior = abs(beta.hat.prior[1:p])+abs(beta.hat.prior[(p+1):(2*p)])
    
    ## Fit the lasso on the data from environment 1, tuning weight of prior
    X_env = cbind(Xs[[1]], X_ks[[1]])
    eval_gamma = function(gamma) {
        penalty = 1*(1-gamma)  + gamma * 1 / (0.05+beta.hat.prior)
        cv_fit = cv.glmnet(X_env, ys[[1]], penalty.factor=rep(penalty,2), dfmax=dfmax, family = family)
        idx.min = which.min(cv_fit$cvlo+cv_fit$cvup)
        err = c(cv_fit$cvlo[idx.min], cv_fit$cvup[idx.min])
        return(err)
    }        
    gamma.seq = seq(0,1,length.out=10) # Sequence of prior weights to consider
    err.seq = sapply(gamma.seq, function(gamma) eval_gamma(gamma))
    
    if(FALSE) {
        tibble(gamma=gamma.seq, low=err.seq[1,], up=err.seq[2,]) %>%
            mutate(mid = (low+up)/2) %>%
            gather(low, up, mid, key="Limit", value="Err") %>%            
            ggplot(aes(x=gamma, y=Err, color=Limit)) +
            geom_point() +
            geom_line()
    }        
    
    ## Find optimal prior weight
    idx.best = which.min(colMeans(err.seq))
    ## Make sure this is significantly better than gamma=0
    if(err.seq[2,idx.best] >= mean(err.seq[,1])) {
        idx.best  = 1
    }
    gamma = gamma.seq[idx.best]
    if(verbose) cat(sprintf("Optimal prior weight: %.3f\n", gamma))
    
    ## Re-fit the lasso on environment 1 with optimally tuned prior
    penalty = 1*(1-gamma)  + gamma * 1 / (0.05+beta.hat.prior)
    cv_fit = cv.glmnet(X_env, ys[[1]], penalty.factor=rep(penalty,2), dfmax=dfmax, family = family)
    
    ## Extract importance measures
    beta.hat = coef(cv_fit, s="lambda.min")[-1][1:(2*p)]
    W1 = abs(beta.hat[1:p]) - abs(beta.hat[(p+1):(2*p)])
    return(W1)
}
