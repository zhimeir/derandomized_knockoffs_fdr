#!/bin/bash
AMP_LIST=(10 11 12 13 14)

for SEED in {1..100}; do
  for AMP in "${AMP_LIST[@]}"; do
    Rscript ../simulations/simulation_highdim.R $SEED $AMP 
  done
done
