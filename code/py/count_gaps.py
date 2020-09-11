""" Count gaps in aligned fasta files """


def main():
    with open(snakemake.input.fasta, "r") as infile:
        n_gaps, total_chars = count_gaps(infile)
    with open(snakemake.output.txt, "w") as outfile:
        outfile.write("n_gaps\ttotal_chars\tdataset\tref_weight\tref_frac\tsample_frac\tseed\n")
        outfile.write(f"{n_gaps}\t{total_chars}\t{snakemake.wildcards.dataset}\t{snakemake.wildcards.ref_weight}\t{snakemake.wildcards.ref_frac}\t{snakemake.wildcards.sample_frac}\t{snakemake.wildcards.seed}\n")


def count_gaps(infile, gap_char="-", header_char=">"):
    seqs_cat = "".join(
        line.strip() for line in infile if not line.startswith(header_char)
    )
    n_gaps = seqs_cat.count(gap_char)
    total_chars = len(seqs_cat)
    return n_gaps, total_chars


main()
