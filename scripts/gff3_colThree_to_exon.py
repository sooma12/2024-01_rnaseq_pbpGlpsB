# gff3_colThree_to_exon.py
"""
Change cds/gene to exon in gff3 genome file for STAR genomeGenerate
gff3 file has the following terms in the 3rd column (index 2):
{'region': 1, 'gene': 3628, 'CDS': 3591, 'rRNA': 18, 'exon': 94, 'tRNA': 72,
'pseudogene': 55, 'sequence_feature': 1, 'riboswitch': 7, 'SRP_RNA': 1,
'ncRNA': 1, 'tmRNA': 1, 'RNase_P_RNA': 1}

"""

GFF_FILE_IN = '/work/geisingerlab/Mark/REF_GENOMES/17978-mff/NZ_CP012004.gff3'
GFF_EXON_FILE_OUT = '/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/ref/NZ_CP012004_col3ex.gff3'

with open(GFF_FILE_IN, 'r') as infile, open(GFF_EXON_FILE_OUT, 'w') as outfile:
    TODO_temp_coldict = {}
    for line in infile:
        split_line = line.rstrip().split(sep='\t')
        # if split_line[0] == 'NZ_CP012004.1':
        #     item = split_line[2]
        #     if item in TODO_temp_coldict:
        #         TODO_temp_coldict[item] += 1
        #     else:
        #         TODO_temp_coldict[item] = 1
        new_line = split_line
        if new_line[0] == 'NZ_CP012004.1':
            new_line[2] = 'exon'
        line_to_write = '\t'.join(new_line)
        outfile.write(f'{line_to_write}\n')
