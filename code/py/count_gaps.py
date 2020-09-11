""" Count gaps in aligned fasta files """


def main():
    with open(snakemake.input.fasta, "r") as infile:
        n_gaps, total_chars = count_gaps(infile)
    with open(snakemake.output.txt, "w") as outfile:
        outfile.write("dataset\tn_gaps\ttotal_chars\n")
        outfile.write(f"{snakemake.wildcards.dataset}\t{n_gaps}\t{total_chars}\n")

def count_gaps(infile, gap_char="-", header_char=">"):
    seqs_cat = ''.join(line.strip() for line in infile if not line.startswith(header_char))
    n_gaps = seqs_cat.count(gap_char)
    total_chars = len(seqs_cat)
    return n_gaps, seqs_cat


main()
