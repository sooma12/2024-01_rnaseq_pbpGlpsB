FASTQ_INDIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/input/fastq/

leftInSuffix="_R1_001.fastq"
rightInSuffix="_R2_001.fastq"
leftOutSuffix="_trimmed_R1_fastq"
rightOutSuffix="_trimmed_R2_fastq"

# Loop through all the left-read fastq files in INDIR
for leftInFile in ${FASTQ_INDIR}*${leftInSuffix}
do
  echo "$leftInFile"
  # Remove path from filename
  pathRemoved="${leftInFile/"$FASTQ_INDIR"/}"
  # Remove left-read suffix, leaving the name of the sample
  sampleName="${pathRemoved/$leftInSuffix/}"
  echo Trimming $sampleName
  # Test with echos; comment this out before final use
  echo "$FASTQ_INDIR$sampleName$leftInSuffix"
  echo "$FASTQ_INDIR$sampleName$rightInSuffix"
  echo "$PAIRED_OUTDIR$sampleName$leftOutSuffix"
  echo "$UNPAIRED_OUTDIR$sampleName$leftOutSuffix"
  echo "$PAIRED_OUTDIR$sampleName$rightOutSuffix"
  echo "$UNPAIRED_OUTDIR$sampleName$rightOutSuffix"
done
