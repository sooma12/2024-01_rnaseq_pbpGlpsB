#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=merge_index_bams
#SBATCH --time=04:00:00
#SBATCH --mem=100G
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/slurm_logs/%x-%j.log
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/slurm_logs/%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=soo.m@northeastern.edu

# Use samtools to merge and index bam files, yielding .merged.bam and .merged.bam.bai files

BAM_IN_DIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data/mapped/
MERGE_DIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data/merged_bams/
BAM_SUFFIX=Aligned.sortedByCoord.out.bam
STRAIN_1=WT
STRAIN_2=DpbpG
STRAIN_3=DlpsB

# Make out directories
mkdir -p MERGE_DIR

# Load modules
echo "Loading anaconda and samtools"
module load anaconda3/2021.11
source activate BINF-12-2021
module load samtools/1.10

# Merge replicates to a single bam file
echo "Merging samples from strain: ${STRAIN_1}"
STRAIN_1_BAM_LIST=${MERGE_DIR}${STRAIN_1}_bams.list
ls ${BAM_IN_DIR}${STRAIN_1}*${BAM_SUFFIX} >${STRAIN_1_BAM_LIST}
samtools merge -b ${STRAIN_1_BAM_LIST} ${STRAIN_1}_merged.bam

echo "Merging samples from strain: ${STRAIN_2}"
STRAIN_2_BAM_LIST=${MERGE_DIR}${STRAIN_2}_bams.list
ls ${BAM_IN_DIR}${STRAIN_2}*${BAM_SUFFIX} >${STRAIN_2_BAM_LIST}
samtools merge -b ${STRAIN_2_BAM_LIST} ${STRAIN_2}_merged.bam

echo "Merging samples from strain: ${STRAIN_3}"
STRAIN_3_BAM_LIST=${MERGE_DIR}${STRAIN_3}_bams.list
ls ${BAM_IN_DIR}${STRAIN_3}*${BAM_SUFFIX} >${STRAIN_3_BAM_LIST}
samtools merge -b ${STRAIN_3_BAM_LIST} ${STRAIN_3}_merged.bam

echo "Indexing merged bam files"
for merged_bam in ${MERGE_DIR}*_merged.bam; do
  echo "Indexing ${merged_bam}"
  samtools index $merged_bam
done
