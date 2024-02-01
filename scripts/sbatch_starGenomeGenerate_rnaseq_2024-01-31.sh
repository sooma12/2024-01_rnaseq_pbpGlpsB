#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=starGenomeGenerate_2024-01-31
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --output=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/slurm_logs/%x-%j.log
#SBATCH --error=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/slurm_logs/%x-%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=soo.m@northeastern.edu

## Usage: sbatch sbatch_starGenomeGenerate_rnaseq_2024-01-31.sh

module load star/2.7.11a

# Variables
GENOME_REF_DIR_OUT=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/ref/
# FASTA_IN=/work/geisingerlab/Mark/REF_GENOMES/17978-mff/NZ_CP012004.fasta
# GTF_IN=/work/geisingerlab/Mark/REF_GENOMES/17978-mff/NZ_CP012004.gff3
NTHREADS=4

# STAR requires the output directory be pre-made
mkdir -p $GENOME_REF_DIR_OUT

# STAR time
# Program recommended `--genomeSAindexNbases 9` after running with default value 14
STAR --runMode genomeGenerate \
--genomeDir $GENOME_REF_DIR_OUT \
--genomeFastaFiles $FASTA_IN \
--sjdbGTFfile $GTF_IN \
--sjdbGTFfeatureExon gene \
--genomeSAindexNbases 9 \
--runThreadN $NTHREADS
