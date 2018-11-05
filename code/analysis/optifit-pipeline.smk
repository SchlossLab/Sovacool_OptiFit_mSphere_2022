import os

configfile: 'code/analysis/config_soil.yaml'
dataset = config['dataset']
input_dir = os.path.join(config['input_dir'], dataset)
output_dir = os.path.join(config['output_dir'], dataset)

# TODO: include snakemake workflows in to download data

def count_input_seqs(infilename):
	with open(infilename, 'r') as infile:
		num_seqs = len([line for line in file if line[0] == '>'])
	return num_seqs
num_seqs = count_input_seqs(os.path.join(input_dir, dataset + '.fasta'))
