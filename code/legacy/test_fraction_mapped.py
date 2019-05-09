#!/usr/local/bin/python3
input_mapped = ['sample.optifit_mcc.list']
input_original = ['sample.count_table']
output = ['sample.fraction_mapped.tsv']
if len(input_mapped) != len(input_original):
    raise ValueError("Unequal number of optifit_mcc.list and count_table files")
with open(output[0], 'w') as output_file:
    output_file.write('original_filename\tmapped_filename\tfraction_mapped\n')
    for mapped_filename, original_filename in zip(input_mapped, input_original):
        with open(original_filename, 'r') as input_file:
            line = next(input_file)  # first column of all lines except first line
            input_samples = set([line.split()[0] for line in input_file])
            print(f"{original_filename}\t{input_samples}")
        with open(mapped_filename, 'r') as mapped_file:
            line = next(mapped_file)
            line = next(mapped_file) # third column onward of second line, each seq in each OTU
            mapped_samples = set(seq for column in line.split()[2:] for seq in column.split(','))
            print(f"{mapped_filename}\t{mapped_samples}")
        fraction_mapped = len(input_samples.intersection(mapped_samples)) / len(input_samples)
        output_file.write(f'{original_filename}\t{mapped_filename}\t{fraction_mapped}\n')
