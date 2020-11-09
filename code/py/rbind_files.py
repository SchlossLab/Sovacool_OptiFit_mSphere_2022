# perform row bind, but don't hold data in memory.
# assumes all tsv files have the same columns.

with open(snakemake.output.tsv, 'w') as outfile:
    header_row = ''
    for infilename in snakemake.input.tsv:
        with open(infilename, 'r') as infile:
            line = next(infile) # throw out header row except first time
            if not header_row:
                header_row = line
                outfile.write(header_row)
            for line in infile:
                outfile.write(line)