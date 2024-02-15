# pipeline_2_trim_align_counts.sh
# Calls scripts to run Trimmomatic, STAR/alignRNA, and subread/featureCounts through SBATCH
# Usage: request a compute node with srun or my ~/request_compute_node.sh script
# FROM WORK_DIR (below): `bash scripts/pipeline_2_trim_align_counts.sh &`
# TODO FOR FINAL USAGE `bash scripts/pipeline_2_trim_align_counts.sh 1>./slurm_logs/last.log 2>./slurm_logs/last.err &`
# gunzipped fastq input data should be in $FASTQ_INDIR.  Don't pass .gz files!

WORK_DIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB
CUR_DATE="$(date '+%Y-%m-%d_%H-%M')"
## TODO Remove "temp" from the paths below.  Just doing this for simplicity.  rm -r ./slurm_logs/temp
LOG_DIR=${WORK_DIR}/slurm_logs/temp/${CUR_DATE}
mkdir -p $LOG_DIR
MAIN_LOG_FILE=${WORK_DIR}/slurm_logs/temp/logfile_${CUR_DATE}.log
touch $MAIN_LOG_FILE

echo "Starting analysis - $(date '+%Y-%m-%d %H:%M:%S')" >>$MAIN_LOG_FILE
echo
echo

echo "Making directory for log files." >>$MAIN_LOG_FILE
mkdir -p $LOG_DIR
echo "Logs will be stored in $LOG_DIR" >>$MAIN_LOG_FILE
echo
echo

echo "Trimming fastq data - $(date '+%Y-%m-%d %H:%M:%S')" >>$MAIN_LOG_FILE
echo "Loading trimmomatic" >>$MAIN_LOG_FILE
module load trimmomatic/0.39
module load oracle_java/jdk1.8.0_181
echo "Running trimmomatic script" >>$MAIN_LOG_FILE
FASTQ_INDIR=$WORK_DIR/input/fastq/
PAIRED_OUTDIR=$WORK_DIR/data/RNA/trimmed/paired/
UNPAIRED_OUTDIR=$WORK_DIR/data/RNA/trimmed/unpaired/
mkdir -p $PAIRED_OUTDIR $UNPAIRED_OUTDIR
## TODO not live version yet.  Echoes.
sbatch --partition=short --job-name=trim_rnaseq --time=02:00:00 -N 1 -n 2 \
--mail-type END --mail-user soo.m@northeastern.edu \
--output=$LOG_DIR/%x-%j.log --error=$LOG_DIR/%x-%j.err \
./scripts/2_align_and_count/trim_quality_length.sh $FASTQ_INDIR $PAIRED_OUTDIR $UNPAIRED_OUTDIR
echo "trimming complete; outputs saved to ${PAIRED_OUTDIR} and ${UNPAIRED_OUTDIR}" >>$MAIN_LOG_FILE
echo
echo


echo "Preparing sample sheet with paired files" >>$MAIN_LOG_FILE
## TODO not live version yet.  Echoes.
SAMPLE_SHEET=${PAIRED_OUTDIR}/sample_sheet.txt
bash ./scripts/2_align_and_count/prep_sample_sheet_for_starAlign.sh $PAIRED_OUTDIR $SAMPLE_SHEET
echo "Sample sheet saved to ${SAMPLE_SHEET}" >>$MAIN_LOG_FILE
echo "Sample sheet contents:" >>$MAIN_LOG_FILE
echo "$(cat $SAMPLE_SHEET)" >>$MAIN_LOG_FILE
echo
echo


## TODO another option for STAR alignment... could run sequentially on each line of SAMPLE_SHEET.
## Do some variable substitution to get the length of SAMPLE SHEET with wc -l <file> | awk '{print $1}'
## Then make an array of numbers from 1-(length)
## sed/awk through the array to grab each line of SAMPLE SHEET

echo "Aligning trimmed RNA reads to genome - $(date '+%Y-%m-%d %H:%M:%S')" >>$MAIN_LOG_FILE
echo "Loading STAR" >>$MAIN_LOG_FILE
module load star/2.7.11a
echo "Running STAR aligner" >>$MAIN_LOG_FILE
GENOME_DIR=$WORK_DIR/ref
ALIGNED_BAM_DIR=$WORK_DIR/data  # Note that the star_align_rna.sh script directs output to ./data/mapped/<dirs>
mkdir -p $ALIGNED_BAM_DIR
## TODO not live version yet.  Echoes.
sbatch --partition short --job-name STARalignRNA --time 04:00:00 \
--array=1-9%10 --ntasks=9 --mem=100G --cpus-per-task=1 \
--mail-type END --mail-user soo.m@northeastern.edu \
--output=$LOG_DIR/%x-%j.log --error=$LOG_DIR/%x-%j.err \
./scripts/star_align_rna.sh $GENOME_DIR $ALIGNED_BAM_DIR $PAIRED_OUTDIR $SAMPLE_SHEET
echo "mapped .bam files saved to ${ALIGNED_BAM_DIR}/mapped/<sample_name>" >>$MAIN_LOG_FILE
echo
echo
exit  ## TODO REMOVE


echo "Generating counts table - $(date '+%Y-%m-%d %H:%M:%S')" >>$MAIN_LOG_FILE
echo "Loading subread package" >>$MAIN_LOG_FILE
module load subread/2.0.6
echo "Running subread.featureCounts" >>$MAIN_LOG_FILE
FEATURECOUNTS_OUT_DIR=$WORK_DIR/data/featurecounts
GTF_REF=$WORK_DIR/ref/NZ_CP012004_transcript2exon.gtf
mkdir -p ${FEATURECOUNTS_OUT_DIR}
## TODO not live version yet.  Echoes.
sbatch --partition short --job-name featureCounts --time 02:00:00 -N 1 -n 4 \
--mail-type END --mail-user soo.m@northeastern.edu \
--output=$LOG_DIR/%x-%j.log --error=$LOG_DIR/%x-%j.err \
./scripts/featurecounts_get_count_table.sh $FEATURECOUNTS_OUT_DIR $GTF_REF $ALIGNED_BAM_DIR
echo "featureCounts output saved to ${FEATURECOUNTS_OUT_DIR}" >>$MAIN_LOG_FILE
echo
echo

echo "Trimmomatic/STAR.alignRNA/featureCount script completed - $(date '+%Y-%m-%d %H:%M:%S')" >>$MAIN_LOG_FILE
