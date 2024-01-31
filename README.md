# 2024-01_rnaseq_pbpGlpsB
Scripts and records for RNA-seq comparisons of WT vs. ∆pbpG vs. ∆lpsB in A. baumannii 17978

## Sample preparation and processing
Followed RNA isolation, reverse transcription, and qPCR protocol to isolate RNA.
Protocol document (followed up to step 22): https://docs.google.com/document/d/1iIDJtlX0jUI4P77D0M6SsOHyaJRPXctn8PhxoLuhKHk/edit
(EG's protocol with MWS's clarifications)
Aliquots of 20 uL of RNA at 100 ng/uL (totaling 2 ug) were made in separate tubes.  These samples were frozen at -80C before shipment to SeqCenter (on ~10 lbs dry ice).

## Library prep and sequencing
Library preparation and sequencing were performed at SeqCenter.  The following is copy/pasted from their methods pdf (provided with data):
```text
Samples were DNAse treated with Invitrogen DNAse (RNAse free). Library preparation was
performed using Illumina’s Stranded Total RNA Prep Ligation with Ribo-Zero Plus kit and 10bp
unique dual indices (UDI). Sequencing was done on a NovaSeq X Plus, producing paired end
150bp reads. Demultiplexing, quality control, and adapter trimming was performed with bcl-
convert (v4.2.4)1. Sequencing statistics are included in the ‘RNA Sequencing Stats.xlsx’ file.
```

## QC
SeqCenter fastq.gz files were verified using md5sum (script `check_md5sums.py`).

gz files were gunzipped.  Resulting fastq files were subjected to QC using FastQC (script `sbatch_fastqc_rnaseq_2024-01-30.sh`).
Versions:
`OpenJDK/19.0.1
fastqc/0.11.9`
Paired fastq files (R1/R2) were also verified to contain the same number of reads using `wc -l *.fastq`

- Overrepresented sequences from FastQC include: ssrA (tRNA), rplB (L12 large subunit ribosomal protein)
- ∆pbpG R1 and R2 files both have overrepresented short sequences ending in N that don't map to AB genome

## Trimming
Adapter trimming was performed by SeqCenter; manual review of FastQC .html files showed no adapter sequences.

Quality trimming and read length filtering were performed using Trimmomatic.


- Adapter trimming done by seqcenter
- Need quality trimming
- Filter out short reads?