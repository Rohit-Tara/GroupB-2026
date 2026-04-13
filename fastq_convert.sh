#!/bin/bash

# Script for converting SRA files to FASTQ formatr using fasterq-dump 
# Using multi-threading and a temporary directory to improve performance 
# Last edited: 01/03/2026. Rohit

# Loads the SRA toolkit using the pathway, so fasterq-dump can be used
export PATH=/home/r/rt334/Documents/Steered_Research_Project/sratoolkit.3.3.0-ubuntu64/bin:$PATH

# Defines the directories
WORKDIR=/scratch/alice/r/rt334/steered_sra # This uses the scratch directory as my main working directory 
OUTDIR=$WORKDIR/fastq_output # Where the FASTQ files will go
TMPDIR=$WORKDIR/tmp # Temporary working direcotry is created to spead things up

# Move to working directory
cd $WORKDIR

# Create output + temp directories
mkdir -p $OUTDIR # Creates the FASTQ output directory
mkdir -p $TMPDIR # Creates the temporary directory for intermediate files

echo "Starting conversion at $(date)"
#Loops thorugh the SRR ID directories
for dir in SRR* # Each SRR folder contains one .sra file 
do
    echo "Processing $dir"
# This converts the SRA into FASTQ
    fasterq-dump "$dir/$dir.sra" \ # fasterq-dump extracts sequencing reads from the .sra file
        --threads 8 \ # Uses 8 CPU threads for faster processing 
        --temp $TMPDIR \ # Uses the created temporary directory to store the intermediate files
        -O $OUTDIR # The output FASTQ files are stored into the fastq_output directory 

done
# The log completion and time
echo "Finished at $(date)"
