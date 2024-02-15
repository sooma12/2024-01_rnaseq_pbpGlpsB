#!/bin/bash
# featurecounts_get_count_table.sh
# Run featurecounts to generate a count table from RNA-seq BAM files (genome-aligned)
# Usage: bash featurecounts_get_count_table.sh <path/to/output> <path/to/reference.gtf> <path/to/STAR/bams>

OUT_DIR=${1}
GTF_REF=${2}
STAR_DIR=${3}

# Run featureCounts on all BAM files from STAR
# TODO REMOVE ECHO!
echo featureCounts \
-a $GTF_REF \
-o $OUT_DIR/counts.txt \
-p \
--countReadPairs \
$STAR_DIR/*.bam
