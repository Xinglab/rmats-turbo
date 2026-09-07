#!/bin/bash

# Path to the first set of BAM files (e.g., sick group)
# Path to second set of BAM files
# Annotation genes in GTF format
# Output directory
# Type of read (paired- or single-end)
# Temporary output directory
# Length of the reads. Required parameter
rmats.py \
  --b1 b1.txt \
  --b2 b2.txt \
  --gtf gtf/Homo_sapiens.Ensembl.GRCh37.75.gtf \
  --od bam_test \
  -t paired \
  --tmp temp_output_rmats \
  --readLength 50 \
  --cstat 0.0001 \
  --libType fr-unstranded