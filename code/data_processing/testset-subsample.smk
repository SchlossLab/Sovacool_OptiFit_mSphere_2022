import os
configfile: 'config_test.yaml'
datasets = config['datasets']
input_dir = config['input_dir']
subsample_size = config['subsample_size']

rule all:
	input:
		expand("{input_dir}/{dataset}/{subsample_size}.{dataset}.{ext}", input_dir=input_dir, subsample_size=subsample_size, dataset=datasets, ext={'fasta', 'accnos', 'count_table'})

rule subset:
	input:
		'{input_dir}/{dataset}/{dataset}.fasta'
	output:
		expand('{{input_dir}}/{{dataset}}/{{subsample_size}}.{{dataset}}.{ext}', ext={'fasta', 'accnos', 'count_table'})
	shell:
		'mothur "#set.dir(output={input_dir}/{dataset}); '
		'sub.sample(inputdir={input_dir}/{dataset}, fasta={dataset}.fasta, size={subsample_size}); '
		'list.seqs(fasta=current); '
		'get.seqs(accnos=current, count={dataset}.count_table); '
		'rename.file(accnos = current, fasta=current, count=current, prefix={subsample_size}.{dataset});"'
