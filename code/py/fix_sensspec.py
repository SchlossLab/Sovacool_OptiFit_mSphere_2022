"""
sensspec files from open-ref optifit have one extra column (`numotus`) than closed-ref.
this is a one-time-only script to fill in `NA` for the `numotus` column in rows where method=closed.
"""
import sys

infilename = sys.argv[1]
outfilename = sys.argv[2]

with open(outfilename, 'w') as outfile:
    with open(infilename, 'r') as infile:
        header = next(infile)
        outfile.write(header)
        for line in infile:
            line_split = line.split('\t')
            if line_split[38] == 'closed':
                line_split.insert(2, 'NA')
                line = '\t'.join(line_split)
            outfile.write(line)
