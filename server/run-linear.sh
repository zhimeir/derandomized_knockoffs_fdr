#!/bin/bash

## parameters
# linear
AMP_LIST=(30 40 50 60 70 80 90)

# Slurm parameters
PART=candes,normal,stat,pilanci     		  # Partition names
MEMO=5G                          # Memory required (10GB)
TIME=00-10:00:00                  # Time required (2h)
CORE=1                          # Cores required (10)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS

# Loop through the configurations
for SEEDA in {1..100}; do
    for AMP in "${AMP_LIST[@]}"; do

        # Script to be run
        SCRIPT="linear.sh $SEEDA $AMP"

        # Define job name for this chromosome
        JOBN="simulation_"$SEEDA"_"$AMP
        OUTF=$LOGS"/"$JOBN".out"
        ERRF=$LOGS"/"$JOBN".err"

        # Assemble slurm order for this job
        ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT

        #Print order
        echo $ORD

        # Submit order
        $ORD
    done
done


