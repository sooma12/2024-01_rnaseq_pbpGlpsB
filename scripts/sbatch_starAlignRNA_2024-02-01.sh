#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=alignRNA_2024-02-01
#SBATCH --time=04:00:00
#SBATCH --array=1-9%10
#SBATCH --ntasks=9
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
# See https://www.biostars.org/p/449164/

# Create .list files with R1 and R2 fastqs.  Sort will put them in same orders, assuming files are paired
find $FASTQ_DIR -maxdepth 1 -name "*.fastq" | grep -e "_R1" | sort > R1.list
find $FASTQ_DIR -maxdepth 1 -name "*.fastq" | grep -e "_R2" | sort > R2.list

if [ -f "${SAMPLE_SHEET}" ] ; then
  rm "${SAMPLE_SHEET}"
fi

# make sample sheet from R1 and R2 files.  Format on each line looks like (space separated):
# WT_1 /path/to/WT_1_R1.fastq /path/to/WT_1_R2.fastq
# from sample sheet, we can access individual items from each line with e.g. `awk '{print $3}' sample_sheet.txt`

paste R1.list R2.list | while read R1 R2 ;
do
    outdir_root=$(echo "${R2}" | cut -f9 -d"/" | cut -f1,2 -d"_")
    sample_line="${outdir_root} ${R1} ${R2}"
    echo "${sample_line}" >> $SAMPLE_SHEET
done

# sed and awk read through the sample sheet and grab each whitespace-separated value
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p $SAMPLE_SHEET |  awk '{print $1}')
r1=$(sed -n "$SLURM_ARRAY_TASK_ID"p $SAMPLE_SHEET |  awk '{print $2}')
r2=$(sed -n "$SLURM_ARRAY_TASK_ID"p $SAMPLE_SHEET |  awk '{print $3}')

echo "Running STAR on files $r1 and $r2"

# STAR time
STAR --runMode alignReads \
--genomeDir $GENOME_DIR \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn $r1 $r2 \
--runThreadN 4 \
--outFileNamePrefix ${OUT_DIR}/mapped/${name}

