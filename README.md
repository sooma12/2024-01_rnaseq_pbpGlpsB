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
First attempt (script `sbatch_trimQualityLength_rnaseq_2024-01-31.sh`) ran Trimmomatic using the following settings:
- PE (paired end data)
- `-phred33` specifies the PHRED quality encoding in the input FASTQ files
- `HEADCROP:0` indicates not to trim any bases from start of reads; unnecessary as reads are not multiplexed and were previously adapter-clipped.
- `LEADING:20` and `TRAILING:20` indicates the minimum quality to keep at the start and end of each read.
- `SLIDINGWINDOW:4:30` specifies the values to use in the sliding window for trimming.  These values will cause Trimmomatic to scan a sliding window of 4 bp and clip when the average quality in that 4-bp window is below 30.
- `MINLEN:36` sets the minimum read length to 36; shorter reads are discarded.

*NOTE: In my first run on 1/31/2024, these trimming parameters resulted in loss of ~15-20% of reads.*
*NOTE2: Can also consider Cutadapt, which was used by TvO lab in one of Eddie's papers*

## Reference Genome

### Genome files
Genome files were downloaded from NCBI Nucleotide using the following command on 1/31/2024:
wget -O <output filename> "<NCBI file URL>"  # Note the "" surrounding the URL!!

The following urls were used to download the fasta and gff3 files for the 17978-mff chromosomes and plasmids pAB1-3:
17978-mff chromosome:
NZ_CP012004.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=NZ_CP012004.1"
NZ_CP012004.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=NZ_CP012004.1"
17978-mff pAB3:
NZ_CP012005.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=NZ_CP012005.1"
NZ_CP012005.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=NZ_CP012005.1"

pAB1 and pAB2 files below were used in EG's EGA83 reference, which I used in my Breseq in Oct 2023:
pAB1:
CP000522.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000522.1"
CP000522.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=CP000522.1"
pAB2:
CP000523.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000523.1"
CP000523.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=CP000523.1"

Files are stored in `/work/geisingerlab/Mark/REF_GENOMES/17978-mff`.

### Formatting genome files for STAR.genomeGenerate
gffread was used to convert the NZ_CP012004.gff3 file to gtf.
https://github.com/gpertea/gffread

gffread installation (to /work/geisingerlab/Mark/software/gffread): 
```bash
cd /work/geisingerlab/Mark/software/gffread
git clone https://github.com/gpertea/gffread
cd gffread
make release
```

This makes a binary file /work/geisingerlab/Mark/software/gffread/gffread.
Then I followed these guides to make gffread a module available to me:
https://researchcomputing.princeton.edu/support/knowledge-base/custom-modules
https://hpc.ncsu.edu/Documents/user_modules.php
(made changes to .bashrc for user-installed modules, and added module files in /work/geisingerlab/Mark/software/modulefiles/gffread)

A gtf file can be generated using: 
`gffread  REF_GENOMES/17978-mff/NZ_CP012004.gff3 -T -o REF_GENOMES/17978-mff/NZ_CP012004.gtf`

The GTF file's third column contains "transcript" and "CDS" information, but STAR relies on the third column containing "exon".
We want the information in "transcript" to replace "exon" for STAR.
Two ways to do this: 
1. Add the flag `--sjdbGTFfeatureExon transcript` to STAR command
2. Or, convert "transcript" to "exon".
I did the latter.  Used custom script `gff3_colThree_to_exon.py` with the following command:
`python gff3_colThree_to_exon.py /work/geisingerlab/Mark/REF_GENOMES/17978-mff/NZ_CP012004.gtf /work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/ref/NZ_CP012004_transcript2exon.gtf transcript`

Finally, STAR's genomeGenerate was run using the commands in `sbatch_starGenomeGenerate_GTF_2024-02-01.sh`

STAR genomeGenerate output is found in `ref/`.

## Alignment
Note, try the next few steps both with quality/length-trimmed data AND with untrimmed data!


