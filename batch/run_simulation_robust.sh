#!/bin/bash
AMP_LIST=(4.5 5 5.5 6 6.5)
N_LIST=(600 700 800 900 1000)
for SEED in {1..100}; do
  for AMP in "${AMP_LIST[@]}"; do
    for N in "${N_LIST[@]}"; do
      Rscript ../simulations/simulation_linear_robust.R $SEED $AMP $N
    done
  done
done
