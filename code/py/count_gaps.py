''' Count gaps in aligned fasta files '''

def main():
    with open(snakemake.input.fasta, 'r') as infile:
        n_gaps = count_gaps(infile)
    with open(snakemake.output.txt, 'w') as outfile:
        outfile.write('dataset\tn_gaps\n')
        outfile.write(f'{snakemake.wildcards.dataset}\t{n_gaps}\n')

def count_gaps(input_filename, gap_char = '-', header_char = '>'):
    return sum(line.count(gap_char) for line in infile if not line.startswith(header_char))

main()