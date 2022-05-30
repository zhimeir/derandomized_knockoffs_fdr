#!/bin/bash
AMP_LIST=(5.5 6 6.5 7 7.5)

for SEED in {1..100}; do
  for AMP in "${AMP_LIST[@]}"; do
    Rscript ../simulations/simulation_multi.R $SEED $AMP 
  done
done
