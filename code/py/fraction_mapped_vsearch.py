#!/usr/local/bin/python3
"""
calculate the fraction of reads that mapped during closed-reference clustering
"""

def count_seqs(infilename):
    with open(infilename, "r") as seq_file:
        num_seqs = sum(1 for line in seq_file if line.startswith('>'))
    return num_seqs

def get_fraction_mapped(num_query_seqs, num_unmapped_seqs):
    print(num_unmapped_seqs, num_query_seqs)
    fraction_mapped = 1 - round(num_unmapped_seqs / num_query_seqs, 3)
    print('frac map:', fraction_mapped)
    assert fraction_mapped <= 1 and fraction_mapped >= 0
    return fraction_mapped

def write_fraction_mapped(wildcards, fraction_mapped):
    method = "closed"#wildcards.method
    ref = "NA" if method == "de_novo" else "gg"
    region = "NA"
    header_line = "dataset\tref\tregion\tmethod\tfraction_mapped\n"
    data_str = f"{wildcards.dataset}\t{ref}\t{region}\t{method}\t{fraction_mapped}\n"
    with open(snakemake.output.txt, "w") as output_file:
        output_file.write(header_line)
        output_file.write(data_str)

def main():
    num_query = count_seqs(snakemake.input.query)
    num_unmapped = count_seqs(snakemake.input.unmapped)
    frac_mapped = get_fraction_mapped(num_query, num_unmapped)
    write_fraction_mapped(snakemake.wildcards, frac_mapped)

main()
