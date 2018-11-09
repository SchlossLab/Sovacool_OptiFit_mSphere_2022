import os
configfile: 'config_test.yaml'
datasets = config['datasets']
input_dir = config['input_dir']
subsample_size = config['subsample_size']

rule all:
	input:
		expand("{input_dir}/{dataset}_{subsample_size}/{dataset}_{subsample_size}.{ext}", input_dir=input_dir, subsample_size=subsample_size, dataset=datasets, ext={'fasta', 'accnos', 'count_table'})

rule subset:
	input:
		'{input_dir}/{dataset}/{dataset}.fasta'
	output:
		expand('{{input_dir}}/{{dataset}}_{{subsample_size}}/{{dataset}}_{{subsample_size}}.{ext}', ext={'fasta', 'accnos', 'count_table'})
	params:
		dataset = '{dataset}'
	shell:
		'mothur "#set.dir(input={input_dir}/{params.dataset}, output={input_dir}/{params.dataset}_{subsample_size}); '
		'sub.sample(fasta={params.dataset}.fasta, size={subsample_size}); '
		'list.seqs(fasta=current); '
		'get.seqs(accnos=current, count={params.dataset}.count_table); '
		'rename.file(accnos = current, fasta=current, count=current, prefix={params.dataset}_{subsample_size});"'
