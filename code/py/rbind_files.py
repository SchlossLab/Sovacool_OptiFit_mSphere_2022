# perform row bind, but don't hold data in memory.
# assumes all tsv files have the same columns and exactly one row.

with open(snakemake.output.tsv, 'w') as outfile:
    for infilename in snakemake.input.tsv:
        with open(infilename, 'r') as infile:
            next(infile) # throw out header row
            outfile.write(next(infile))