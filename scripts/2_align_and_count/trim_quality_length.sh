#!/bin/bash
# trim_quality_length.sh
# Run trimmomatic on RNA-seq fastq files
# MWS Feb. 14th, 2024
# Usage: bash trim_quality_length.sh <path/to/fastq/inputs> <path/to/output/pairedfiles> <path/to/output/unpairedfiles>

# adapter trimming is done by seqcenter; the following is unnecessary but kept for future use.
# path to NU Discovery cluster's Trimmomatic program folder with Illumina adapters
# PATH_TO_TRIMMOMATIC="/shared/centos7/anaconda3/2021.11/envs/BINF-12-2021/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2"

## Initialize variables to contain file suffixes and output paths
leftInSuffix="_R1_001.fastq"
rightInSuffix="_R2_001.fastq"
leftOutSuffix="_trimmed_R1.fastq"
rightOutSuffix="_trimmed_R2.fastq"

FASTQ_INDIR=${1}
PAIRED_OUTDIR=${2}
UNPAIRED_OUTDIR=${3}

# set variables for trimmomatic
HEADCROP=0
LEADING=20
TRAILING=20
SLIDINGWINDOW='4:15'
MINLEN=36

echo 'trimmomatic parameters:'
echo "HEADCROP: $HEADCROP"
echo "LEADING: $LEADING"
echo "TRAILING: $TRAILING"
echo "SLIDINGWINDOW: $SLIDINGWINDOW"
echo "MINLEN: $MINLEN"

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
  # TODO: REMOVE ECHO BEFORE RUNNING FOR REAL
  echo java -jar /shared/centos7/trimmomatic/0.39/trimmomatic-0.39.jar PE \
  -threads 1 -phred33 \
  $FASTQ_INDIR$sampleName$leftInSuffix \
  $FASTQ_INDIR$sampleName$rightInSuffix \
  $PAIRED_OUTDIR$sampleName$leftOutSuffix \
  $UNPAIRED_OUTDIR$sampleName$leftOutSuffix \
  $PAIRED_OUTDIR$sampleName$rightOutSuffix \
  $UNPAIRED_OUTDIR$sampleName$rightOutSuffix \
  HEADCROP:${HEADCROP} \
  LEADING:${LEADING} \
  TRAILING:${TRAILING} \
  SLIDINGWINDOW:${SLIDINGWINDOW} \
  MINLEN:${MINLEN}
  # ILLUMINACLIP:$PATH_TO_TRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10
done
