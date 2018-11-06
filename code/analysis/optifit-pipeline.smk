""" Benchmarking the OptiFit algorithm """

import math
import os

configfile: 'config_soil.yaml'

input_dir = os.path.join(config['input_dir'], config['dataset'])
output_dir = os.path.join(config['output_dir'], config['dataset'])
dataset = config['dataset'] if not config['subsample_test'] else '.'.join([str(config['subsample_size']), config['dataset']])

def count_input_seqs(infilename):
	with open(infilename, 'r') as infile:
		num_seqs = len([line for line in infile if line[0] == '>'])
	return num_seqs
num_seqs = count_input_seqs(os.path.join(input_dir, dataset + '.fasta'))
weight = config['weight']
sizes = [math.floor(num_seqs * i/100) for i in range(5,10,5)]
iters = range(config['iterations'])
reps = range(config['replicates'])

rule all:
	input:
		expand('{input_dir}/{dataset}.{ext}', dataset=dataset, input_dir=input_dir, ext={'fasta', 'dist', 'count_table'}),
		expand('{output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/sample.opti_mcc.{ext}', output_dir=output_dir, dataset=dataset, weight=weight, size=sizes, iter=iters, rep=reps, method={'open', 'closed'}, printref={'t', 'f'}, ext={'list', 'steps', 'sensspec'})

rule get_dists:
	input:
		'{input_dir}/{{dataset}}.fasta'.format(input_dir=input_dir)
	output:
		'{input_dir}/{{dataset}}.dist'.format(input_dir=input_dir)
	shell:
		'mothur "#set.dir(input={input_dir}, output={input_dir}); '
		'dist.seqs(fasta={dataset}.fasta, cutoff=0.03);"'

rule split_weighted_subsample:
	input:
		count="{input_dir}/{{dataset}}.count_table".format(input_dir=input_dir),
		dist="{input_dir}/{{dataset}}.dist".format(input_dir=input_dir)
	params:
		size="{size}",
		weight="{weight}"
	output:
		"{output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.accnos"
	script:
		"weighted_subsample.R"

rule prep_weighted_subsample:
	input:
		fasta="{input_dir}/{{dataset}}.fasta".format(input_dir=input_dir),
		count="{input_dir}/{{dataset}}.count_table".format(input_dir=input_dir),
		dist="{input_dir}/{{dataset}}.dist".format(input_dir=input_dir),
		accnos="{output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.accnos"
	output:
		fasta="{output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.fasta",
		count="{output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.count_table",
		dist="{output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.dist"
	shell:
		'mothur "#set.seed(seed={iter}}); '
		'set.dir(output={output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/); '
		'get.seqs(accnos={input.accnos}, fasta={input.fasta}, count={input.count_table}); '
		'get.dists(column={input.dist}, accnos=current); '
		'rename.file(fasta=current, count=current, accnos = current, column=current, prefix=sample)"'

rule prep_reference_from_dataset:
	input:
		fasta="{input_dir}/{dataset}.fasta",
		count="{input_dir}/{dataset}.count_table",
		dist="{input_dir}/{dataset}.dist",
		accnos="{output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.accnos"
	output:
		expand("{{output_dir}}/dataset-as-reference/{{dataset}}_weight-{{weight}}_size-{{size}}_i-{{iter}}/r-{{rep}}/reference.{ext}", ext={"accnos", 'count_table', 'dist', 'fasta'})
	shell:
		'mothur "#set.seed(seed={iter}); '
		'set.dir(output={output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/); '
		'remove.seqs(fasta={input.fasta}, count={input.count}, accnos={input.accnos}); '
		'list.seqs(fasta=current); '
		'get.dists(column={input.dist}, accnos=current); '
		'rename.file(fasta=current, count=current, column=current, accnos=current, prefix=reference)"'

rule cluster_samples:
	input:
		count="{output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.count_table",
		column="{output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.dist"
	output:
		expand('{{output_dir}}/dataset-as-reference/{{dataset}}_weight-{{weight}}_size-{{size}}_i-{{iter}}/r-{{rep}}/sample.opti_mcc.{ext}', ext={'list', 'steps', 'sensspec'})
	shell:
		'mothur "#set.seed(seed={rep}); '
		'set.dir(output={output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/); '
		'cluster(column={input.column}, count={input.count})"'

rule cluster_reference:
	input:
		count="{output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/reference.count_table",
		column="{output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/reference.dist"
	output:
		expand('{{output_dir}}/dataset-as-reference/{{dataset}}_weight-{{weight}}_size-{{size}}_i-{{iter}}/r-{{rep}}/reference.opti_mcc.{ext}', ext={'list', 'steps', 'sensspec'})
	shell:
		'mothur "#set.seed(seed={rep}); '
		'set.dir(output={output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/); '
		'cluster(column={input.column}, count={input.count})"'

rule fit:
	input:
		reflist='reference.opti_mcc.list',
		refcolumn='reference.dist',
		refcount='reference.count_table',
		reffasta='reference.fasta',
		fasta='sample.fasta',
		count='sample.count_table',
		column='sample.dist',
	output:
		expand('{{output_dir}}/dataset-as-reference/{{dataset}}_weight-{{weight}}_size-{{size}}_i-{{iter}}/r-{{rep}}/method-{{method}}_printref-{{printref}}/sample.opti_mcc.{ext}', ext={'list', 'steps', 'sensspec'})
	shell:
		'mothur "set.seed(seed={rep}); '
		'set.dir(input={output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/, output={output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/); '
		'cluster.fit(reflist={input.reflist}, refcolumn={input.refcolumn}, refcount={input.refcount}, reffasta={input.reffasta}, fasta={input.fasta}, count={input.count}, column={input.column}, printref={printref}, method={method}); '
