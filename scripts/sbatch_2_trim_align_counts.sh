#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=trim_align_counts_2024-02-14
#SBATCH --time=08:00:00
#SBATCH --array=1-9%10
#SBATCH --ntasks=9
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=soo.m@northeastern.edu

# usage:  `sbatch sbatch_2_trim_align_counts.sh
# gunzipped fastq input data should be in $FASTQ_INDIR.  Don't pass .gz files!

WD=/Users/mws/Documents/geisinger_lab_research/bioinformatics_in_acinetobacter/rnaSeq/2024-01_rnaseq_pbpGlpsB
LOG_DIR=${WD}/slurm_logs/$SLURM_JOB_ID

echo "Starting analysis - $(date)"

echo "Setting working directory: $(WD)"
cd $WD || exit
echo "Making directory for log files."
mkdir -p $LOG_DIR
echo "Logs will be stored in $LOG_DIR"

echo "Trimming fastq data - $(date)"
echo "Loading trimmomatic"
module load trimmomatic/0.39
module load oracle_java/jdk1.8.0_181
echo "Running trimmomatic script"
FASTQ_INDIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/input/fastq/
PAIRED_OUTDIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data/RNA/trimmed/paired/
UNPAIRED_OUTDIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data/RNA/trimmed/unpaired/
mkdir -p $PAIRED_OUTDIR $UNPAIRED_OUTDIR
## TODO not live version yet.  Echoes.
bash ./scripts/2_align_and_count/trim_quality_length.sh $FASTQ_INDIR $PAIRED_OUTDIR $UNPAIRED_OUTDIR \
1>$LOG_DIR/$SLURM_JOB_NAME-$SLURM_JOB_ID-trim.log 2>$LOG_DIR/$SLURM_JOB_NAME-$SLURM_JOB_ID-trim.err
echo "trimming complete; outputs saved to ${PAIRED_OUTDIR} and ${UNPAIRED_OUTDIR}"

echo "Preparing sample sheet with paired files"
## TODO not live version yet.  Echoes.
SAMPLE_SHEET=${PAIRED_OUTDIR}/sample_sheet.txt
bash ./scripts/2_align_and_count/prep_sample_sheet_for_starAlign.sh $PAIRED_OUTDIR $SAMPLE_SHEET
echo "Sample sheet saved to ${SAMPLE_SHEET}"

echo "Aligning trimmed RNA reads to genome - $(date)"
echo "Loading STAR"
module load star/2.7.11a
echo "Running STAR aligner"
GENOME_DIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/ref
ALIGNED_BAM_DIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data
mkdir -p $ALIGNED_BAM_DIR
## TODO not live version yet.  Echoes.
bash ./scripts/star_align_rna.sh $GENOME_DIR $ALIGNED_BAM_DIR $PAIRED_OUTDIR $SAMPLE_SHEET \
1>$LOG_DIR/$SLURM_JOB_NAME-$SLURM_JOB_ID-align.log 2>$LOG_DIR/$SLURM_JOB_NAME-$SLURM_JOB_ID-align.err
echo "mapped .bam files saved to ${ALIGNED_BAM_DIR}/mapped/<sample_name>"

echo "Generating counts table - $(date)"
echo "Loading subread package"
module load subread/2.0.6
echo "Running subread.featureCounts"
FEATURECOUNTS_OUT_DIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data/featurecounts
GTF_REF=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/ref/NZ_CP012004_transcript2exon.gtf
mkdir -p ${FEATURECOUNTS_OUT_DIR}
## TODO not live version yet.  Echoes.
bash ./scripts/featurecounts_get_count_table.sh $FEATURECOUNTS_OUT_DIR $GTF_REF $ALIGNED_BAM_DIR \
1>$LOG_DIR/$SLURM_JOB_NAME-$SLURM_JOB_ID-count.log 2>$LOG_DIR/$SLURM_JOB_NAME-$SLURM_JOB_ID-count.err
echo "featureCounts output saved to ${FEATURECOUNTS_OUT_DIR}"

echo "Trimmomatic/STAR.alignRNA/featureCount script completed - $(date)"
