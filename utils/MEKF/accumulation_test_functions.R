##############################################################################
## These functions implement 
##   the Accumulation Test methods from the paper:
## Ang Li & Rina Foygel Barber,
##   "Accumulation tests for FDR control
##      in ordered hypothesis testing"
## Available from http://arxiv.org/abs/1505.07352
## (Several methods from other papers also implemented,
##        as noted below - see citations in paper)
##############################################################################

##############################################################################
## HingeExp method,
##    i.e. an accumulation test with the HingeExp function:
##       h(p) = C * log(1/(C*(1-p))) * (p>1-1/C)
##############################################################################
create_HingeExp_function = function(C=2){
	function(x){C*log(1/(C*(1-x)))*(x>1-1/C)}
}
HingeExp = function(pvals,alpha=0.2,C=2,output_type='khat'){
	AccumulationTest(pvals,create_HingeExp_function(C),alpha=alpha,output_type=output_type,check_integrate_to_one=FALSE)
}
##############################################################################


##############################################################################
## ForwardStop method (G'Sell et al 2013),
##    i.e. an accumulation test with the ForwardStop function:
##       h(p) = log(1/(1-p))
##############################################################################
create_ForwardStop_function = function(){
	function(x){log(1/(1-x))}
}
ForwardStop = function(pvals,alpha=0.2,output_type='khat'){
	AccumulationTest(pvals,create_ForwardStop_function(),alpha=alpha,output_type=output_type,check_integrate_to_one=FALSE)
}
##############################################################################


##############################################################################
## SeqStep method (Barber&Candes 2015),
##    i.e. an accumulation test with the step function:
##       h(p) = C * (p>1-1/C)
##############################################################################
create_SeqStep_function = function(C=2){
	function(x){C*(x>1-1/C)}
}
SeqStep = function(pvals,alpha=0.2,C=2,output_type='khat'){
	AccumulationTest(pvals,create_SeqStep_function(C),alpha=alpha,output_type=output_type,check_integrate_to_one=FALSE)
}
####################################################################################


##############################################################################
## SeqStep+ method (Barber&Candes 2015),
##    i.e. an accumulation test with the step function:
##       h(p) = C * (p>1-1/C)
##         & with the conservative correction
##              for estimating FDR
##############################################################################
SeqStepPlus = function(pvals,alpha=0.2,C=2,output_type='khat'){
	AccumulationTest(pvals,create_SeqStep_function(C),alpha=alpha,numerator_plus=C,denominator_plus=1,output_type=output_type,check_integrate_to_one=FALSE)
}
##############################################################################


##############################################################################
## Accumulation test for a generic function "hfun"
##############################################################################
AccumulationTest = function(pvals,hfun,alpha=0.2,numerator_plus=0,denominator_plus=0,output_type='khat',check_integrate_to_one=TRUE){
	
	# check for valid arguments
	check_inputs = CheckInputs_AccumulationTest(pvals,hfun,alpha,numerator_plus,denominator_plus,output_type,check_integrate_to_one)
	if(length(check_inputs)>0){
		stop(check_inputs)
	}
	
	# perform the test
	n=length(pvals)
	FDPest=(numerator_plus+cumsum(unlist(lapply(pvals,hfun))))/(denominator_plus+1:n)
	FDPest_vs_alpha=(FDPest%*%t(rep(1,length(alpha)))<=rep(1,n)%*%t(alpha))
	findlast=function(x){max(c(0,which(x)))}
	khat=apply(FDPest_vs_alpha,2,findlast)
	if(output_type=='khat'){		
		return(khat)
	}else{if(output_type=='FDPest'){
		return(FDPest)
	}else{
		output=list()
		output$FDPest=FDPest
		output$khat=khat
		return(output)
	}}
}
##############################################################################



##############################################################################
## Check inputs for AccumulationTest
##############################################################################
CheckInputs_AccumulationTest = function(pvals,hfun,alpha=0.2,numerator_plus=0,denominator_plus=0,output_type='khat', check_integrate_to_one){
	# check_integrate_to_one should be logical
	if(!is.logical(check_integrate_to_one)){
		return('check_integrate_to_one must be logical')
	}
	# check that pvals and alpha are each sequences with values in [0,1]
	if(!is.numeric(pvals) || !is.vector(pvals) || min(pvals)<0 || max(pvals)>1){
		return('pvals must be a number or numeric vector with values in [0,1]')
	}
	n=length(pvals)
	
	if(!is.numeric(alpha) || !is.vector(alpha) || min(alpha)<0 || max(alpha)>1){
		return('alpha must be a number or numeric vector with values in [0,1]')
	}	
	
	# check that hfun is a function that gives a nonnegative value for each pvalue
	if(!is.function(hfun)){
		return('hfun must be a function')
	}
	if(!is.numeric(try(hfun(pvals),silent=TRUE)) || any(is.na(try(hfun(pvals),silent=TRUE)))){
		return('The function hfun must take as input a vector of p-values in [0,1], and return a vector of nonnegative numbers of the same length')
	}
	if(length(hfun(pvals))!=length(pvals)){
		return('The function hfun must take as input a vector of p-values in [0,1], and return a vector of nonnegative numbers of the same length')
	}
	if(any(hfun(pvals)<0)){
		return('The function hfun must take as input a vector of p-values in [0,1], and return a vector of nonnegative numbers of the same length')
	}
	if(check_integrate_to_one){
		if(abs(integrate(hfun,0,1)$value-1)>1e-2){
			return('The function hfun must have expectation 1 when applied to a uniform variable p~Uniform[0,1] (set check_integrate_to_one=FALSE to override)')
		}
	}
	
	# check that numerator_plus and denominator_plus are numbers
	if(!is.numeric(numerator_plus) || !is.numeric(denominator_plus) || length(numerator_plus)!=1 || length(denominator_plus)!=1){
		return('numerator_plus and denominator_plus must each be a scalar')
	}
	
	# check that output_type is in {'khat', 'FDPest', 'both'}
	if(!is.element(output_type,c('khat','FDPest','both'))){
		return('Invalid output type: choose from "khat", "FDPest", or "both"')
	}
	
	return(NULL)
}
##############################################################################

