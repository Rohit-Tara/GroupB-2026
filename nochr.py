# Script to remove (chr) prefix from chromosome names in the GTF file 
# Last edited on: 16/03/2026
# Converts GENCODE annotation to match genome index format as it has no chr 

# The code removes (chr) from the first column while keeping header lines the same
awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} {sub(/^chr/, "", $1); print}' \
gencode.v22.annotation.gtf > gencode.v22.nochr.annotation.gtf
