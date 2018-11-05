import os
configfile: 'config_soil.yaml'
dataset = config['dataset']
input_dir = os.path.join(config['input_dir'], dataset)

def count_input_seqs(infilename):
	with open(infilename, 'r') as infile:
		num_seqs = len([line for line in infile if line[0] == '>'])
	return num_seqs
num_seqs = count_input_seqs(os.path.join(input_dir, dataset + '.fasta'))

rule subset:
	input:
		os.path.join(input_dir, dataset + '.fasta')
	output:
	shell:
		
