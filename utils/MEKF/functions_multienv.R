################################################################
## These functions run the multi-environment knockoff analysis
################################################################
library(knockoff)

## Report discoveries from testing consistent conditional associations
## Inputs: 
##    W_matrix: a matrix of multi-environment knockoff statistics
##    q: FDR threshold
## Outputs:
##    Variables selected by the multi-environment knockoff filter
invariant_model = function(W_matrix, q = 0.1){
    W_sign = 2*((apply(W_matrix, 1, min) > 0) - 0.5)
    W_prod = abs(apply(W_matrix, 1, prod))
    W_eff = W_sign * W_prod
    thres = knockoff.threshold(W_eff, fdr=q)
    selected = which(W_eff >= thres)
    return(selected)
}

## Find the sign of entries in a vector
## Inputs: 
##    vec: a vector
## Outputs:
##    The sign of the each entry: -1 for negative, +1 for positive, and +1/-1
##              with probability 0.5 for 0. 
sign_fun = function(vec){
    l = length(vec)
    return((rbinom(l,1,0.5)*(vec == 0) + (vec > 0))*2 - 1)
}

## Compute partial conjunction multi-environment p-values
## Inputs: 
##    vec: a vector of knockoff statistics
##    r: the number of nonnull environments required
##    randomness: randomness = 0 computes the p-values as in (16) in the paper without the Uj term
##                randomness = 1 computes the p-values as in (19) in the paper without the Uj term
##                randomness = 2 computes the p-values as in (16) in the paper with the Uj term
##                randomness = 3 computes the p-values as in (19) in the paper with the Uj term      
## Outputs:
##    Partial conjunction multi-environment p-values
p_value = function(vec, r, randomness){
    l = length(vec)
    minus = sum(vec < 0)
    num_zero = sum(vec == 0)
    
    if(randomness == 3){
        minus = minus + rbinom(1,num_zero,0.5)
        num_zero = 0
    }
    
    a0 = l - r + 1
    a = l - r + 1 - num_zero
    if (a > 0){
        discrete_p = pbinom(minus,a,0.5)
        discrete_p_l = pbinom(minus - 1,a,0.5)
    } else {
        discrete_p = 1
        discrete_p_l = 0
    }
    p_cont = runif(1,discrete_p_l,discrete_p)
    steps = pbinom(0:a0,a0,0.5)
    p_disc = steps[min(which(steps >= p_cont))]
    if(randomness == 0){
        return (discrete_p)
    }else if(randomness == 1){
        return (p_disc)
    }else if(randomness >= 2){
        return (p_cont)
    }
}

## Combine multi-environment knockoff statistics
## Inputs: 
##    vec: a vector of knockoff statistics
##    r: the number of nonnull environments required
## Outputs:
##    invariant statistics
combine_mag = function(vec, r){
    vec_abs = abs(vec)
    vec_sort = -sort(-vec_abs, partial = r)[1:r]
    return(prod(vec_sort))
}


## Report discoveries from testing partically consistent conditional associations
## Inputs: 
##    W_matrix: a matrix of multi-environment knockoff statistics
##    r: the number of nonnull environments required
##    q: FDR threshold
##    method: variable section method; "seqstep" stands for Selective Seqstep+, "accumulation" stands for Accumulation test
##    c: threshold c in Selective Seqstep+
##    randomness: randomness = 0 computes the p-values as in (16) in the paper without the Uj term
##                randomness = 1 computes the p-values as in (19) in the paper without the Uj term
##                randomness = 2 computes the p-values as in (16) in the paper with the Uj term
##                randomness = 3 computes the p-values as in (19) in the paper with the Uj term      
## Outputs:
##    Variables selected by the multi-environment knockoff filter
partial_conjunction = function(W_matrix, r, q = 0.1, method = "seqstep", c = 0.6, randomness = 2){
    pvals = apply(W_matrix, 1, p_value, r = r, randomness = randomness)
    W_sign = 2*(pvals <= c) - 1
    W_mag = apply(W_matrix, 1, combine_mag, r = r)
    if (method == "seqstep"){
        q_tilde = q*(1-c)/c
        W_eff = W_sign * W_mag
        thres = knockoff.threshold(W_eff, fdr=q_tilde)
        selected = which(W_eff >= thres)
        return (selected)
    } else if(method == "accumulation"){
        hfun = create_HingeExp_function(C=2)
        pvals_sorted = pvals[order(-W_mag)]
        num_selecteda = AccumulationTest(pvals_sorted , hfun, alpha=q)
        selecteda = NULL
        if (num_selecteda > 0){
            selecteda = order(-W_mag)[1:num_selecteda]
        }
        return(sort(selecteda))
    }
    return (NULL)
}