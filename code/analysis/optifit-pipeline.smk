""" Benchmarking the OptiFit algorithm """

import math
import os

configfile: 'config_test.yaml'

input_dir = config['input_dir']
output_dir = config['output_dir']
weight = config['weight']
iters = range(config['iterations'])
reps = range(config['replicates'])
methods = {'open', 'closed'}
printrefs = {'t', 'f'}

class Dataset:
	def __init__(self, name, subsample_size, fasta_filename):
		self.name = name
		self.subsample = subsample_size
		self.num_seqs = self.count_input_seqs(fasta_filename)
		self.sizes = [math.floor(self.num_seqs * i/100) for i in range(5,11,5)]

	def count_input_seqs(self, fasta_filename):
		with open(fasta_filename, 'r') as infile:
			num_seqs = len([line for line in infile if line[0] == '>'])
		return num_seqs

datasets = dict()
for dataset_name in config['datasets']:
	subsample_size = config['subsample_size'] if config['subsample_test'] else ''
	fasta_filename = os.path.join(input_dir, dataset_name, '{dataset}{subsample}.fasta'.format(dataset=dataset_name, subsample=subsample_size))
	datasets[name] = Dataset(dataset_name, subsample_size, fasta_filename)

rule all:
	input:
		['{input_dir}/{dataset}/{dataset}{subsample}.{ext}'.format( dataset=datasets[name], input_dir=input_dir, subsample=datasets[name].subsample, ext=extension) for name in datasets for extension in {'fasta', 'dist', 'count_table'}],
		["{output_dir}/dataset-as-reference/{dataset}/{dataset}{subsample}_weight-{weight}_size-{size}_i-{iter}/sample.{ext}".format(output_dir=output_dir, subsample=datasets[name].subsample, dataset=name, weight=weight, size=size, iter=iter, ext=ext) for name in datasets for size in datasets[name].sizes for iter in iters for ext in {'fasta','count_table', 'dist'}],
		['{output_dir}/dataset-as-reference/{dataset}/{dataset}{subsample}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/sample.opti_mcc.{ext}'.format( output_dir=output_dir, subsample=datasets[name].subsample, dataset=name, weight=weight, size=size, iter=iter, rep=rep, method=method, printref=printref, ext=ext) for name in datasets for size in datasets[name].sizes for iter in iters for rep in reps for method in methods for printref in printrefs for ext in {'list', 'steps', 'sensspec'}]

rule get_dists:
	input:
		'{input_dir}/{dataset}/{dataset}{subsample}.fasta'
	output:
		'{input_dir}/{dataset}/{dataset}{subsample}.dist'
	shell:
		'mothur "#set.dir(input={input_dir}/{dataset}/, output={input_dir}/{dataset}/); '
		'dist.seqs(fasta={dataset}{subsample}.fasta, cutoff=0.03);"'

rule split_weighted_subsample:
	input:
		count="{input_dir}/{dataset}{subsample}.count_table",
		dist="{input_dir}/{dataset}{subsample}.dist"
	params:
		size="{size}",
		weight="{weight}"
	output:
		"{output_dir}/dataset-as-reference/{dataset}/{dataset}{subsample}_weight-{weight}_size-{size}_i-{iter}/sample.accnos"
	script:
		"weighted_subsample.R"

rule prep_weighted_subsample:
	input:
		fasta="{input_dir}/{dataset}{subsample}.fasta",
		count="{input_dir}/{dataset}{subsample}.count_table",
		dist="{input_dir}/{dataset}{subsample}.dist",
		accnos="{output_dir}/dataset-as-reference/{dataset}/{dataset}{subsample}_weight-{weight}_size-{size}_i-{iter}/sample.accnos"
	output:
		expand("{{output_dir}}/dataset-as-reference/{{dataset}}/{{dataset}}{{subsample}}_weight-{{weight}}_size-{{size}}_i-{{iter}}/sample.{ext}", ext=['fasta','count_table', 'dist'])
	shell:
		'mothur "#set.seed(seed={iter}}); '
		'set.dir(output={output_dir}/dataset-as-reference/{dataset}{subsample}_weight-{weight}_size-{size}_i-{iter}/); '
		'get.seqs(accnos={input.accnos}, fasta={input.fasta}, count={input.count_table}); '
		'get.dists(column={input.dist}, accnos=current); '
		'rename.file(fasta=current, count=current, accnos = current, column=current, prefix=sample)"'

rule prep_reference_from_dataset:
	input:
		fasta="{input_dir}/{dataset}{subsample}.fasta",
		count="{input_dir}/{dataset}{subsample}.count_table",
		dist="{input_dir}/{dataset}{subsample}.dist",
		accnos="{output_dir}/dataset-as-reference/{dataset}/{dataset}{subsample}_weight-{weight}_size-{size}_i-{iter}/sample.accnos"
	output:
		expand("{{output_dir}}/dataset-as-reference/{{dataset}}/{{dataset}}{{subsample}}_weight-{{weight}}_size-{{size}}_i-{{iter}}/r-{{rep}}/reference.{ext}", ext={"accnos", 'count_table', 'dist', 'fasta'})
	shell:
		'mothur "#set.seed(seed={iter}); '
		'set.dir(output={output_dir}/dataset-as-reference/{dataset}/{dataset}{subsample}_weight-{weight}_size-{size}_i-{iter}/); '
		'remove.seqs(fasta={input.fasta}, count={input.count}, accnos={input.accnos}); '
		'list.seqs(fasta=current); '
		'get.dists(column={input.dist}, accnos=current); '
		'rename.file(fasta=current, count=current, column=current, accnos=current, prefix=reference)"'

rule cluster_samples:
	input:
		count="{output_dir}/dataset-as-reference/{dataset}/{dataset}{subsample}_weight-{weight}_size-{size}_i-{iter}/sample.count_table",
		column="{output_dir}/dataset-as-reference/{dataset}/{dataset}{subsample}_weight-{weight}_size-{size}_i-{iter}/sample.dist"
	output:
		expand('{{output_dir}}/dataset-as-reference/{{dataset}}/{{dataset}}{{subsample}}_weight-{{weight}}_size-{{size}}_i-{{iter}}/r-{{rep}}/sample.opti_mcc.{ext}', ext={'list', 'steps', 'sensspec'})
	shell:
		'mothur "#set.seed(seed={rep}); '
		'set.dir(output={output_dir}/dataset-as-reference/{dataset}/{dataset}{subsample}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/); '
		'cluster(column={input.column}, count={input.count})"'

rule cluster_reference:
	input:
		count="{output_dir}/dataset-as-reference/{dataset}/{dataset}{subsample}_weight-{weight}_size-{size}_i-{iter}/reference.count_table",
		column="{output_dir}/dataset-as-reference/{dataset}/{dataset}{subsample}_weight-{weight}_size-{size}_i-{iter}/reference.dist"
	output:
		expand('{{output_dir}}/dataset-as-reference/{{dataset}}/{{dataset}}{{subsample}}_weight-{{weight}}_size-{{size}}_i-{{iter}}/r-{{rep}}/reference.opti_mcc.{ext}', ext={'list', 'steps', 'sensspec'})
	shell:
		'mothur "#set.seed(seed={rep}); '
		'set.dir(output={output_dir}/dataset-as-reference/{dataset}/{dataset}{subsample}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/); '
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
		expand('{{output_dir}}/dataset-as-reference/{{dataset}}/{{dataset}}{{subsample}}_weight-{{weight}}_size-{{size}}_i-{{iter}}/r-{{rep}}/method-{{method}}_printref-{{printref}}/sample.opti_mcc.{ext}', ext={'list', 'steps', 'sensspec'})
	shell:
		'mothur "set.seed(seed={rep}); '
		'set.dir(input={output_dir}/dataset-as-reference/{dataset}/{dataset}{subsample}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/, output={output_dir}/dataset-as-reference/{dataset}/{dataset}{subsample}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/); '
		'cluster.fit(reflist={input.reflist}, refcolumn={input.refcolumn}, refcount={input.refcount}, reffasta={input.reffasta}, fasta={input.fasta}, count={input.count}, column={input.column}, printref={printref}, method={method}); '
