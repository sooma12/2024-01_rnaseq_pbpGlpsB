#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=trimmo_rnaseq_2024-01-30
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/slurm_logs/%x-%j.log
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/slurm_logs/%x-%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=soo.m@northeastern.edu

## Usage: sbatch sbatch_trimQualityLength_rnaseq_2024-01-31.sh
## Outputs paired and unpaired trimmed reads to data/RNA/trimmed/paired and data/RNA/trimmed/unpaired

# Load trimmomatic
module load trimmomatic/0.39

# Directories for inputs and outputs.  NOTE: MUST END WITH / FOR SHELL SUBSTITUTION BELOW TO WORK
FASTQ_INDIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/input/fastq/
PAIRED_OUTDIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data/RNA/trimmed/paired/
UNPAIRED_OUTDIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data/RNA/trimmed/unpaired/

# path to NU Discovery cluster's Trimmomatic program folder with Illumina adapters
# PATH_TO_TRIMMOMATIC="/shared/centos7/anaconda3/2021.11/envs/BINF-12-2021/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2"

## Initialize variables to contain file suffixes and output paths
leftInSuffix="_R1_001.fastq"
rightInSuffix="_R2_001.fastq"
leftOutSuffix="_trimmed_R1_fastq"
rightOutSuffix="_trimmed_R2_fastq"

mkdir -p $PAIRED_OUTDIR $UNPAIRED_OUTDIR

# Loop through all the left-read fastq files in INDIR
for leftInFile in ${FASTQ_INDIR}*${leftInSuffix}
do
  # Remove path from filename
  pathRemoved="${leftInFile/"$FASTQ_INDIR"/}"
  # Remove left-read suffix, leaving the name of the sample
  sampleName="${pathRemoved/$leftInSuffix/}"
  echo Trimming $sampleName
  # Test with echos; comment this out before final use
#  echo "$FASTQ_INDIR$sampleName$leftInSuffix"
#  echo "$FASTQ_INDIR$sampleName$rightInSuffix"
#  echo "$PAIRED_OUTDIR$sampleName$leftOutSuffix"
#  echo "$UNPAIRED_OUTDIR$sampleName$leftOutSuffix"
#  echo "$PAIRED_OUTDIR$sampleName$rightOutSuffix"
#  echo "$UNPAIRED_OUTDIR$sampleName$rightOutSuffix"
  # Use sample name derived from shell replacement to trim left AND right reads
  trimmomatic PE \
  -threads 1 -phred33 \
  $FASTQ_INDIR$sampleName$leftInSuffix \
  $FASTQ_INDIR$sampleName$rightInSuffix \
  $PAIRED_OUTDIR$sampleName$leftOutSuffix \
  $UNPAIRED_OUTDIR$sampleName$leftOutSuffix \
  $PAIRED_OUTDIR$sampleName$rightOutSuffix \
  $UNPAIRED_OUTDIR$sampleName$rightOutSuffix \
  HEADCROP:0 \
  LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:36 \
  # ILLUMINACLIP:$PATH_TO_TRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 \
done
