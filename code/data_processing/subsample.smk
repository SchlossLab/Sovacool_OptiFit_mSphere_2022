import os
configfile: 'config_soil.yaml'
dataset = config['dataset']
trimsize = config['trimsize']
input_dir = os.path.join(config['input_dir'], dataset)
'''
def count_input_seqs(infilename):
	with open(infilename, 'r') as infile:
		num_seqs = len([line for line in infile if line[0] == '>'])
	return num_seqs
num_seqs = count_input_seqs(os.path.join(input_dir, dataset + '.fasta'))
'''
rule all:
	input:
		expand("{input_dir}/subsets/{trimsize}.{dataset}.{ext}", input_dir=input_dir, trimsize=trimsize, dataset=dataset, ext={'fasta', 'accnos', 'count_table'})

rule subset:
	input:
		'{input_dir}/{dataset}.fasta'
	output:
		expand('{{input_dir}}/subsets/{{trimsize}}.{{dataset}}.{ext}', ext={'fasta', 'accnos', 'count_table'})
	shell:
		'mothur "#set.dir(output={input_dir}/subsets/); '
		'sub.sample(inputdir={input_dir}, fasta={dataset}.fasta, size={trimsize}); '
		'list.seqs(fasta=current); '
		'get.seqs(accnos=current, count={dataset}.count_table); '
		'rename.file(accnos = current, fasta=current, count=current, prefix={trimsize}.{dataset});"'
