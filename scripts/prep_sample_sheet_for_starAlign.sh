# prep_sample_sheet_for_starAlign.sh
# Makes a sample_sheet.txt containing sample ID and R1 and R2 filepaths
# Example output line:  WT_1 /path/to/WT_1_R1.fastq /path/to/WT_1_R2.fastq

## IMPORTANT: CHANGE FASTQ_DIR IN BOTH THIS SCRIPT AND IN ALIGNRNA SBATCH SCRIPT!

FASTQ_DIR=/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/input/fastq
SAMPLE_SHEET=${FASTQ_DIR}/sample_sheet.txt

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

rm R1.list R2.list
