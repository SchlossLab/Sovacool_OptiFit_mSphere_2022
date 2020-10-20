from collections import OrderedDict

input_filename = 'feature-table.tsv'
output_filename = 'otu_table.list'

with open(input_filename, 'r') as infile:
    next(infile) # header line
    next(infile) # column names
    otus = OrderedDict({line.split()[0].strip(): line.split()[1].strip()
                        for line in infile})

with open(output_filename, 'w') as outfile:
    header = f"label\tnumOtus\t{otu + '\t' for otu in otus.keys()}"
    outfile.write(header)

