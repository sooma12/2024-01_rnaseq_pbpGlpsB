#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=subreadFeatureCounts_2024-02-02
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/slurm_logs/%x-%j.log
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/slurm_logs/%x-%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=soo.m@northeastern.edu

# Usage: sbatch sbatch_featurecounts_2024-02-02.sh

OUT_DIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data/featurecounts
GTF_REF=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/ref/NZ_CP012004_transcript2exon.gtf
STAR_DIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data/mapped

echo "Loading environment and tools"
module load anaconda3/2021.05
eval "$(conda shell.bash hook)"
conda activate /work/geisingerlab/conda_env/subread

mkdir -p $OUT_DIR

# Run featureCounts on all BAM files from STAR
featureCounts \
-a $GTF_REF \
-f GTF \
-o $OUT_DIR \
-p \
--countReadPairs \
$STAR_DIR/*.bam
