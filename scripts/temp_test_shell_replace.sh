FASTQ_INDIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/input/fastq

echo $FASTQ_INDIR
leftInSuffix="_R1_001.fastq"

# Loop through all the left-read fastq files in INDIR
for leftInFile in ${FASTQ_INDIR}*${leftInSuffix}
do
  echo "$leftInFile"
done
