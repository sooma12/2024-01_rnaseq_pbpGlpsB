#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=alignRNA_2024-02-01
#SBATCH --time=04:00:00
#SBATCH --array=1-9%10
#SBATCH --ntasks=9
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/slurm_logs/alignRNA.%N.%j.log
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/slurm_logs/alignRNA.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=soo.m@northeastern.edu

## Usage: sbatch sbatch_starAlignRNA_2024-02-01.sh
## Updated 2/1/2024

module load star/2.7.11a

# Variables
GENOME_DIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/ref
OUT_DIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data
# fastq directory... non-quality-trimmed reads are in input/
# quality-trimmed reads are in data/RNA/trimmed/paired and data/RNA/trimmed/unpaired
FASTQ_DIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/input/fastq
SAMPLE_SHEET=${FASTQ_DIR}/sample_sheet.txt

# STAR requires the output directory be pre-made
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
STAR --runMode alignReads \
--genomeDir $GENOME_DIR \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn $r1 $r2 \
--runThreadN 4 \
--alignIntronMax 1 \
--limitBAMsortRAM 5000000000 \
--outFileNamePrefix ${OUT_DIR}/mapped/${name}

