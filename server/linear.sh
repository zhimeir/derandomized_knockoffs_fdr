#!/bin/bash

SEEDA=$1
AMP=$2
ml gmp
ml mpfr
ml gcc
module load R/3.5

Rscript ../simulation/simulation_linear.R $SEEDA $SEEDB $AMP

