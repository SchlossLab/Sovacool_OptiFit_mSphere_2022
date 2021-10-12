#!/usr/local/bin/python3
"""
calculate the fraction of reads that mapped during closed-reference clustering
"""
import numpy as np

def count_seqs(infilename):
    with open(infilename, "r") as seq_file:
        num_seqs = sum(1 for line in seq_file if line.startswith('>'))
    return num_seqs

def get_fraction_mapped(num_query_seqs, num_unmapped_seqs):
    fraction_mapped = 1 - round(num_unmapped_seqs / num_query_seqs, 3)
    assert fraction_mapped <= 1 and fraction_mapped >= 0
    return fraction_mapped


def main():
    num_query = count_seqs(snakemake.input.query)
    num_unmapped = count_seqs(snakemake.input.unmapped)
    frac_mapped_fasta = get_fraction_mapped(num_query, num_unmapped)

    with open(snakemake.input.tsv, 'r') as tsv_file:
        header_line = next(tsv_file).split('\t')
        print(header_line)
        assert header_line[4].strip() == 'fraction_mapped'
        data_line = next(tsv_file).split('\t')
        frac_mapped_list = float(data_line[4].strip())

    print('frac_mapped_fasta: ', frac_mapped_fasta)
    print('frac_mapped_list:  ', frac_mapped_list)

    is_close = np.isclose(frac_mapped_fasta, frac_mapped_list)
    print('is_close: ', is_close)
    #assert is_close

    with open(snakemake.output.txt, 'w') as outfile:
        outfile.write('frac_mapped_fasta\tfrac_mapped_list\tis_close\n')
        outfile.write(f"{frac_mapped_fasta}\t{frac_mapped_list}\t{is_close}\n")

main()
