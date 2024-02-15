#!/bin/bash
# star_align_rna.sh
# Run STAR alignReads to align FASTQ RNA sequences to reference genome
# Usage: bash star_align_rna.sh <path/to/reference/genome> <path/to/output/dir> <path/to/sample/sheet>

# Variables
GENOME_DIR=${1}
OUT_DIR=${2}
SAMPLE_SHEET=${3}

mkdir -p $OUT_DIR

# Array job!  Used my sample sheet technique from 2023-11-13 breseq for this.
# NOTE!!! sample sheet prep is moved to prep_sample_sheet_for_starAlign.sh.
# Run that script before this one!

# sed and awk read through the sample sheet and grab each whitespace-separated value
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p $SAMPLE_SHEET |  awk '{print $1}')
r1=$(sed -n "$SLURM_ARRAY_TASK_ID"p $SAMPLE_SHEET |  awk '{print $2}')
r2=$(sed -n "$SLURM_ARRAY_TASK_ID"p $SAMPLE_SHEET |  awk '{print $3}')

echo "Running STAR on files $r1 and $r2"

# STAR time
# note alignIntronMax=1 to disallow introns for bacteria
# BAM sorting seems to be memory intensive.  Set an available amount of memory using limitBAMsortRAM
# TODO REMOVE ECHO!
STAR --runMode alignReads \
--genomeDir $GENOME_DIR \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn $r1 $r2 \
--runThreadN 4 \
--alignIntronMax 1 \
--limitBAMsortRAM 5000000000 \
--outFileNamePrefix ${OUT_DIR}/mapped/${name}

