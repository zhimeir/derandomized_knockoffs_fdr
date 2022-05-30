#!/bin/bash
AMP_LIST=(16 19 22 25 28)

for SEED in {1..100}; do
  for AMP in "${AMP_LIST[@]}"; do
    Rscript ../simulations/simulation_binom.R $SEED $AMP 
  done
done
