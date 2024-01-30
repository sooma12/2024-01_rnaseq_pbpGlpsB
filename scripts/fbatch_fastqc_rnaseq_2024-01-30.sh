#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=fastqc_rnaseq_2024-01-30
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/slurm_logs/%x-%j.output
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/slurm_logs/%x-%j.error
#SBATCH --mail-type=END
#SBATCH --mail-user=soo.m@northeastern.edu

echo "Starting fastqc SBATCH script $(date)"

echo "Loading environment and tools"
#fastqc requires OpenJDK/19.0.1
module load OpenJDK/19.0.1
module load fastqc/0.11.9

FASTQDIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/input
OUT_DIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/fastqc_output
SCRIPT_DIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/scripts

mkdir -p $FASTQDIR $OUT_DIR $SCRIPT_DIR

echo "Running fastqc in directory $FASTQDIR"
fastqc $FASTQDIR/*.fastq

echo "Cleaning up logs and output files"
mkdir -p $SCRIPT_DIR/logs
mv $SCRIPT_DIR/fastqc_rnaseq* $SCRIPT_DIR/logs
mkdir -p $OUT_DIR/fastqc_html $OUT_DIR/fastqc_zip
mv $FASTQDIR/*fastqc.html $OUT_DIR/fastqc_html
mv $FASTQDIR/*fastqc.zip $OUT_DIR/fastqc_zip



