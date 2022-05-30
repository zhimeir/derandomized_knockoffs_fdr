#!/bin/bash
AMP_LIST=(4.5 5 5.5 6 6.5)
V_LIST=(1 6 11 16 21)
for SEED in {1..100}; do
  for AMP in "${AMP_LIST[@]}"; do
    for V in "${V_LIST[@]}"; do
      Rscript ../simulations/simulation_linear_robust.R $SEED $AMP $V
    done
  done
done
