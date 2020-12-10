#!/usr/local/bin/python3
"""
calculate the fraction of reads that mapped during closed-reference clustering
"""


def parse_seqs(infilename):
    with open(infilename, "r") as count_file:
        if not infilename.endswith(".accnos"):
            next(count_file)  # drop header line if it's not an accnos file
        # first column of all lines
        seqs = {line.split("\t")[0].strip() for line in count_file}
    return seqs


def main():
    query_seqs = parse_seqs(snakemake.input.query)
    ref_seqs = parse_seqs(snakemake.input.ref)
    mapped_seqs = parse_seqs(snakemake.input.mapped)

    print("map==ref", mapped_seqs == ref_seqs)
    print("map==query", mapped_seqs == query_seqs)

    print("len(ref)", len(ref_seqs))
    print("len(query)", len(query_seqs))
    print("len(map)", len(mapped_seqs))

    print("ref int query", len(query_seqs.intersection(ref_seqs)))
    print("map int query", len(mapped_seqs.intersection(query_seqs)))
    print("map int ref", len(mapped_seqs.intersection(ref_seqs)))

    # remove reference seqs from mapped seqs
    mapped_seqs = mapped_seqs - ref_seqs
    print("len(map) after remove ref seqs", len(mapped_seqs))

    # mapped reads should not contain reads not in query
    assert not (mapped_seqs - query_seqs)

    fraction_mapped = round(
        len(mapped_seqs.intersection(query_seqs)) / len(query_seqs), 3
    )
    assert fraction_mapped <= 1 and fraction_mapped >= 0

    wildcards = snakemake.wildcards
    method = wildcards.method
    ref = "NA" if method == "de_novo" else "gg"
    region = "NA"
    header_line = "dataset\tref\tregion\tmethod\tfraction_mapped\n"
    data_str = f"{wildcards.dataset}\t{ref}\t{region}\t{method}\t{fraction_mapped}\n"
    with open(snakemake.output.txt, "w") as output_file:
        output_file.write(header_line)
        output_file.write(data_str)


main()
