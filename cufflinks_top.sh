#!/bin/bash

# Script for running Cufflinks on all TopHat BAM outputs 
# Last edited: 27/03/2026 

# SLURM jobs settings
#SBATCH --job-name=cufflinks_all
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/cufflinks_%j.out
#SBATCH --error=logs/cufflinks_%j.err

# Defining all the paths
BASE="/scratch/alice/r/rt334/steered_sra" #This will be used as a base directory
BAMDIR="${BASE}/tophat_output" # Directory containing TopHat BAM outputs
OUTDIR="${BASE}/cufflinks_output_v2" # An output direcrtory that will contain the cufflinks results. 

# Loading Cufflinks on ALICE HPC
module load cufflinks

#Creates the output directories
mkdir -p "$OUTDIR" logs # Ensures output and the log directories exist with this code

# Created a for loop that loops through all the BAM files, each BAM file being one sample that has been processed by TopHat
for BAM in ${BAMDIR}/SRR*_tophat/accepted_hits.bam
do
    CELL=$(basename $(dirname "$BAM") _tophat) # This extract samples names (SRR ID) from the file path

    echo "Processing $CELL"
# This now runks Cufflinks
    cufflinks \
        -p 8 \ # Uses 8 CPU threads
        -G ${BASE}/gencode.v22.nochr.annotation.gtf \ # Reference gene annotation that being the (nochr that I fixed)
        -o ${OUTDIR}/${CELL}_cufflinks \ # Output for this sample
        "$BAM" # Input BAM file which are the aligned reads

done
# The final message letting me know that all the cufflink outputs have been created
echo "ALL DONE"
