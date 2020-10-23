#!/usr/local/bin/python3
"""
calculate the fraction of reads that mapped during closed-reference clustering
"""


def main():  # TODO: write to log file
    with open(snakemake.input.count, "r") as count_file:
        if not snakemake.input.count.endswith(".accnos"):
            next(count_file)  # drop header line if it's not an accnos file
        # first column of all lines
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
    dsrr = wildcards.dataset_ref_region.split("_")
    ref = dsrr[1] if len(dsrr) > 1 else "NA"
    region = "_".join(dsrr[2:]) if len(dsrr) > 1 else "NA"
    header_line = "dataset\tref\tregion\tmethod\tfraction_mapped\n"
    data_str = (
        f"{wildcards.dataset}\t{ref}\t{region}\t{wildcards.method}\t{fraction_mapped}\n"
    )
    with open(snakemake.output.txt, "w") as output_file:
        output_file.write(header_line)
        output_file.write(data_str)


main()
