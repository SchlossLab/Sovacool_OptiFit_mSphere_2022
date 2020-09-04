#!/usr/local/bin/python3
"""
calculate the fraction of reads that mapped during closed-reference clustering
"""


def main():  # TODO: write to log file
    with open(snakemake.input.count, "r") as count_file:
        # first column of all lines except first line
        if not snakemake.input.count.endswith(".accnos"):
            next(count_file)
        all_reads = {line.split("\t")[0].strip() for line in count_file}
    with open(snakemake.input.list, "r") as list_file:
        # third column onward of second line, each seq in each OTU delimited by comma
        next(list_file)
        line = next(list_file)
        mapped_reads = {seq for column in line.split()[2:] for seq in column.split(",")}
    fraction_mapped = round(
        len(mapped_reads.intersection(all_reads)) / len(all_reads), 3
    )
    wildcards = snakemake.wildcards
    with open(snakemake.output.txt, "w") as output_file:
        data_str = (
            f"{wildcards.dataset}\t{wildcards.seed}\t{wildcards.method}\t{wildcards.printref}\t{fraction_mapped}\t{wildcards.sample_frac}\t{wildcards.ref_frac}\t{wildcards.ref_weight}\n"
            if "sample_frac" in wildcards.__dict__.keys()
            else f"{wildcards.dataset}\t{wildcards.ref}\t{wildcards.region}\t{wildcards.seed}\t{wildcards.method}\t{wildcards.printref}\t{fraction_mapped}\n"
        )
        output_file.write(data_str)


main()
