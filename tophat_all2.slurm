#!/bin/bash

# Script for running TopHat on all FASTQ samples with the limited parallel jobs 
# Last edited: 29/03/2026. Rohit
# Uses paired-end FASTQ files to generate BAM alignments, running up to 4 jobs at once

#SBATCH --job-name=tophat_all2
#SBATCH --partition=medium
#SBATCH --cpus-per-task=32   # 4 jobs × 8 threads each
#SBATCH --mem=64G            # increase memory accordingly
#SBATCH --time=48:00:00
#SBATCH --output=logs/tophat_%j.out
#SBATCH --error=logs/tophat_%j.err

# Defining the main directories
BASE="/scratch/alice/r/rt334/steered_sra"
FASTQ="/scratch/alice/r/rt334/steered_sra/fastq_output"
OUTDIR="/scratch/alice/r/rt334/steered_sra/tophat_output"
# Path to TopHat executable
TOPHAT="/home/r/rt334/miniforge3/envs/tophat_env/bin/tophat"

# Creates output and log directories
mkdir -p "$OUTDIR" logs

# Moves to base directories 
cd "$BASE"

# Set maximum number of parallel jobs
MAX_JOBS=4
RUNNING=0
# Loop through paired-end FASTQ files
for R1 in ${FASTQ}/*_1.fastq
do
    CELL=$(basename "$R1" _1.fastq) # Extracts the sample name
    R2="${FASTQ}/${CELL}_2.fastq"
    # skip samples if already completed
    if [[ -f "${OUTDIR}/${CELL}_tophat/accepted_hits.bam" ]]; then
        echo "Skipping $CELL"
        continue
    fi
    echo "Processing $CELL"
    (
        $TOPHAT \ # Runs TopHat alignment 
            -p 8 \
            -G ${BASE}/gencode.v22.nochr.annotation.gtf \
            -o ${OUTDIR}/${CELL}_tophat \
            ${BASE}/GRCh38_index \
            "$R1" \
            "$R2"
        if [[ -f "${OUTDIR}/${CELL}_tophat/accepted_hits.bam" ]]; then # Checks if alignment is completed succesfully
            echo "SUCCESS: $CELL"
        else
            echo "ERROR: $CELL failed"
            exit 1
        fi
    ) &
    ((RUNNING++))
    # limit to 4 concurrent jobs
    if [[ $RUNNING -ge $MAX_JOBS ]]; then
        wait -n
        ((RUNNING--))
    fi
done
# wait for remaining jobs
wait
echo "ALL DONE" # The final message to confirm completion 
