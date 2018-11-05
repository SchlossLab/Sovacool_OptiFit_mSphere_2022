import os
configfile: 'config_soil.yaml'
dataset = config['dataset']
subsample_size = config['subsample_size']
input_dir = os.path.join(config['input_dir'], dataset)

rule all:
	input:
		expand("{input_dir}/{subsample_size}.{dataset}.{ext}", input_dir=input_dir, subsample_size=subsample_size, dataset=dataset, ext={'fasta', 'accnos', 'count_table'})

rule subset:
	input:
		'{input_dir}/{dataset}.fasta'
	output:
		expand('{{input_dir}}/{{subsample_size}}.{{dataset}}.{ext}', ext={'fasta', 'accnos', 'count_table'})
	shell:
		'mothur "#set.dir(output={input_dir}); '
		'sub.sample(inputdir={input_dir}, fasta={dataset}.fasta, size={subsample_size}); '
		'list.seqs(fasta=current); '
		'get.seqs(accnos=current, count={dataset}.count_table); '
		'rename.file(accnos = current, fasta=current, count=current, prefix={subsample_size}.{dataset});"'
