#!/bin/bash

# Script to retry the failed SRA to FASTQ conversions 
# Last edited: 01/03/2026. Rohit 
# Reads SRR IDs from missing.txt and reprocesses only failed samples 

# Load SRA toolkit so fasterq-dump and validation can be used
export PATH=/home/r/rt334/Documents/Steered_Research_Project/sratoolkit.3.3.0-ubuntu64/bin:$PATH

# Defines the directories 
WORKDIR=/scratch/alice/r/rt334/steered_sra
OUTDIR=$WORKDIR/fastq_output
TMPDIR=$WORKDIR/tmp_retry

# Moves to the working directory 
cd $WORKDIR

# Creating output, temporary and log directories 
mkdir -p $OUTDIR
mkdir -p $TMPDIR
mkdir -p $WORKDIR/retry_logs

# Log start and show missing samples 
echo "Retry started: $(date)"
echo "Missing list:"
cat missing.txt
# Loops through the failed SRR IDs
while read -r srr; do
  echo "----"
  echo "Retrying $srr at $(date)"

  # Remove partial FASTQ outputs from failed runs
  rm -f "$OUTDIR/${srr}_1.fastq" "$OUTDIR/${srr}_2.fastq"

  # Validate .sra before converting
  vdb-validate "$srr/$srr.sra" || { echo "VALIDATION_FAILED $srr"; continue; }

  # Convert SRA to FASTQ (paired-end) 
  fasterq-dump "$srr/$srr.sra" \
    --threads 8 \
    --split-files \
    --temp "$TMPDIR" \
    -O "$OUTDIR" \
    > "$WORKDIR/retry_logs/${srr}.out" 2> "$WORKDIR/retry_logs/${srr}.err"

  # Confirms that both FASTQ files were created and not empty
  if [[ -s "$OUTDIR/${srr}_1.fastq" && -s "$OUTDIR/${srr}_2.fastq" ]]; then
    echo "OK $srr"
  else
    echo "FAILED_OUTPUT $srr"
  fi

done < missing.txt
# Completed the rerun and generated the FASTQs at the specific time it has given 
echo "Retry finished: $(date)"
