#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=srna_featurecounts
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/slurm_logs/%x-%j.log
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/slurm_logs/%x-%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=soo.m@northeastern.edu

# Usage: sbatch sbatch_srna_featurecounts_2024-02-28.sh

OUT_DIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data/featurecounts
GTF_REF=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/ref/17978-mff_sRNAs_2024-02-28.gtf
STAR_DIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data/mapped

echo "Loading environment and tools"
module load subread/2.0.6

mkdir -p $OUT_DIR

# Run featureCounts on all BAM files from STAR
featureCounts \
-a $GTF_REF \
-o $OUT_DIR/srna_counts.txt \
-p \
-t sRNA \
--countReadPairs \
$STAR_DIR/*.bam
