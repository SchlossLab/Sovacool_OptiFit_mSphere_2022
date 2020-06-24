#!/usr/local/bin/python3
"""
calculate the fraction of reads that mapped during closed-reference clustering
"""

def main():
    with open(snakemake.input.count, 'r') as count_file:
        # first column of all lines except first line
        next(count_file)
        all_reads = {line.split('\t')[0] for line in count_file}
    with open(snakemake.input.list, 'r') as list_file:
        # third column onward of second line, each seq in each OTU delimited by comma
        next(list_file)
        line = next(list_file)
        mapped_reads = {seq for column in line.split()[2:] for seq in column.split(',')}
    fraction_mapped = round(len(mapped_reads.intersection(all_reads)) / len(all_reads), 3)
    wildcards = snakemake.wildcards
    with open(snakemake.output.txt, 'w') as output_file:
        output_file.write(f"{wildcards.dataset}\t{wildcards.ref}\t{wildcards.region}\t{wildcards.seed}\t{wildcards.method}\t{wildcards.printref}\t{fraction_mapped}\n")


main()