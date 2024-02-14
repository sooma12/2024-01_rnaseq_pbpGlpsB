# get_a1s_locus_tags_from_acx.py
"""
Given a csv file with a single column containing ACX60 locus tags (from 17978 genome),
use a correspondence table to write out a new csv file with A1S locus tags
For KEGG database!

Based on the following script from Sept 2023:
# make_genelist_with_operons.py
# MWS 9/22/2023
"""

import csv


root_path = "/Users/mws/Documents/geisinger_lab_research/bioinformatics_in_acinetobacter/rnaSeq/2024-01_rnaseq_pbpGlpsB/"
annotations_csv = root_path + "data/labels_TnSeqEssentiality.csv"
infile = root_path + "data/pbpG_results_2024-02-07.csv"
outfile = root_path + "data/pbpG_a1s_degs.csv"


def get_acx_to_a1s_dict(annotations_csv_path):
    """ Read csv file containing annotations; return dictionary of ACX_locus: A1S_locus pairs"""

    acx_to_a1s_dict = {}
    with open(annotations_csv_path, 'r', encoding='utf-8-sig') as annos:  # Had weird characters at start of lines; encoding issue.
        anno_reader = csv.DictReader(annos)
        for row in anno_reader:
            acx_to_a1s_dict[row['RefSeq_locus']] = row['17978_A1S_locus']

    return acx_to_a1s_dict


def acxlist_to_a1s(acx_list: list, acx_a1s_dict: dict):
    """ Take a list of genes with ACX60_ locus tags; return a list converted to A1S locus tags """
    a1s_list = []
    for locus in acx_list:
        try:
            a1s_list.append(acx_a1s_dict[locus])
        except KeyError:  #  skip genes w/o a1s match
            continue
    return a1s_list


def write_genelist_to_csv(out_file_path, regulated_gene_list, column_name):
    """ Take a filepath and a list of genes; write a CSV file just containing one column with list of genes """

    with open(out_file_path, 'w', newline='') as out_fh:
        fieldnames = [f'{column_name}']
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames)
        writer.writeheader()
        for gene in regulated_gene_list:
            if len(gene) > 0:
                writer.writerow({f'{column_name}': gene})


if __name__ == "__main__":
    # Get dictionary to convert ACX locus tags to A1S locus tags
    acx_a1s_dict = get_acx_to_a1s_dict(annotations_csv_path=annotations_csv)

    # Get a list of target genes
    targets_list_acx = []
    with open(infile, 'r') as in_fh:
        line = next(in_fh)
        for line in in_fh:
            line = line.strip()
            split_line = line.split(sep=',')
            locus = split_line[0]
            locus = locus.replace("gene-", "")
            locus = locus.replace('"', '')
            # Keep log2FC > 2 and padj < 0.1
            if split_line[2] == 'NA' or split_line[5] == 'NA':
                continue
            if abs(float(split_line[2])) > 2 and float(split_line[5]) < 0.1:
                targets_list_acx.append(locus)

    # Convert gene list to a1s locus tags
    target_a1s_list = acxlist_to_a1s(targets_list_acx, acx_a1s_dict=acx_a1s_dict)

    # Write gene lists to csvs
    for locus in target_a1s_list:
        if len(locus) < 1:
            target_a1s_list.remove(locus)
    write_genelist_to_csv(outfile, target_a1s_list, column_name='a1s_locus')
