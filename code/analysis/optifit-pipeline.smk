""" Benchmarking the OptiFit algorithm """

import math
import os

configfile: 'config_test.yaml'

input_dir = config['input_dir']
output_dir = config['output_dir']
weights = config['weights']
iters = range(config['iterations'])
reps = range(config['replicates'])
methods = {'open', 'closed'}
printrefs = {'t', 'f'}

class Dataset:
	def __init__(self, name, fasta_filename):
		self.name = name
		self.num_seqs = self.count_input_seqs(fasta_filename)
		self.sizes = {i: math.floor(self.num_seqs * i/100) for i in range(50, 56, 5)}  # { % seqs: num seqs in cut }

	def count_input_seqs(self, fasta_filename):
		with open(fasta_filename, 'r') as infile:
			num_seqs = len([line for line in infile if line[0] == '>'])
		return num_seqs

datasets = dict()
for dataset_name in config['datasets']:
	dataset_name = "{}_{}".format(dataset_name, config['subsample_size']) if config['subsample_test'] else dataset_name
	fasta_filename = os.path.join(input_dir, dataset_name, '{}.fasta'.format(dataset_name))
	datasets[dataset_name] = Dataset(dataset_name, fasta_filename)

rule all:
	input:
		['{input_dir}/{dataset}/{dataset}.{ext}'.format( dataset=name, input_dir=input_dir, ext=extension) for name in datasets for extension in {'fasta', 'dist', 'count_table'}],
		["{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.{ext}".format(output_dir=output_dir, dataset=name, weight=weights, size=size, iter=iter, ext=ext) for name in datasets for size in datasets[name].sizes for iter in iters for ext in {'fasta','count_table', 'dist'}],
		['{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/sample.opti_mcc.{ext}'.format( output_dir=output_dir, dataset=name, weight=weights, size=size, iter=iter, rep=rep, method=method, printref=printref, ext=ext) for name in datasets for size in datasets[name].sizes for iter in iters for rep in reps for method in methods for printref in printrefs for ext in {'list', 'steps', 'sensspec'}]

rule get_dists:
	input:
		'{input_dir}/{dataset}/{dataset}.fasta'
	output:
		'{input_dir}/{dataset}/{dataset}.dist'
	params:
		output_dir='{input_dir}/{dataset}/{dataset}/'
	shell:
		'mothur "#set.dir(output={params.output_dir}); '
		'dist.seqs(fasta={input[0]}, cutoff=0.03);"'

rule split_weighted_subsample:
	input:
		count="{input_dir}/{{dataset}}/{{dataset}}.count_table".format(input_dir=input_dir),
		dist="{input_dir}/{{dataset}}/{{dataset}}.dist".format(input_dir=input_dir)
	params:
		num_seqs=datasets['{dataset}'].sizes["{size}"]
		weight="{weight}"
	output:
		"{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.accnos"
	script:
		"weighted_subsample.R"

rule prep_weighted_subsample:
	input:
		fasta="{input_dir}/{{dataset}}/{{dataset}}.fasta".format(input_dir=input_dir),
		count="{input_dir}/{{dataset}}/{{dataset}}.count_table".format(input_dir=input_dir),
		dist="{input_dir}/{{dataset}}/{{dataset}}.dist".format(input_dir=input_dir),
		accnos="{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.accnos"
	output:
		expand("{{output_dir}}/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{{weight}}_size-{{size}}_i-{{iter}}/sample.{ext}", ext=['fasta','count_table', 'dist'])
	params:
		output_dir="{output_dir}/dataset-as-reference/{dataset}_weight-{weight}_size-{size}_i-{iter}/",
		iter="{iter}"
	shell:
		'mothur "#set.seed(seed={params.iter}); '
		'set.dir(output={params.output_dir}); '
		'get.seqs(accnos={input.accnos}, fasta={input.fasta}, count={input.count}); '
		'get.dists(column={input.dist}, accnos=current); '
		'rename.file(fasta=current, count=current, accnos = current, column=current, prefix=sample)"'

rule prep_reference_from_dataset:
	input:
		fasta="{input_dir}/{{dataset}}/{{dataset}}.fasta".format(input_dir=input_dir),
		count="{input_dir}/{{dataset}}/{{dataset}}.count_table".format(input_dir=input_dir),
		dist="{input_dir}/{{dataset}}/{{dataset}}.dist".format(input_dir=input_dir),
		accnos="{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.accnos"
	output:
		expand("{{output_dir}}/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{{weight}}_size-{{size}}_i-{{iter}}/reference.{ext}", ext={"accnos", 'count_table', 'dist', 'fasta'})
	wildcard_constraints:
		size="\d+",
		iter="\d+"
	params:
		output_dir="{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/",
		iter="{iter}"
	shell:
		'mothur "#set.seed(seed={params.iter}); '
		'set.dir(output={params.output_dir}); '
		'remove.seqs(fasta={input.fasta}, count={input.count}, accnos={input.accnos}); '
		'list.seqs(fasta=current); '
		'get.dists(column={input.dist}, accnos=current); '
		'rename.file(fasta=current, count=current, column=current, accnos=current, prefix=reference)"'

rule cluster_samples:
	input:
		count="{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.count_table",
		column="{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.dist"
	output:
		expand('{{output_dir}}/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{{weight}}_size-{{size}}_i-{{iter}}/r-{{rep}}/sample.opti_mcc.{ext}', ext={'list', 'steps', 'sensspec'})
	wildcard_constraints:
		size="\d+",
		iter="\d+",
		rep="\d+"
	params:
		output_dir='{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/',
		rep="{rep}"
	shell:
		'mothur "#set.seed(seed={params.rep}); '
		'set.dir(output={params.output_dir}); '
		'cluster(column={input.column}, count={input.count})"'

rule cluster_reference:
	input:
		count="{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/reference.count_table",
		column="{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/reference.dist"
	output:
		expand('{{output_dir}}/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{{weight}}_size-{{size}}_i-{{iter}}/r-{{rep}}/reference.opti_mcc.{ext}', ext={'list', 'steps', 'sensspec'})
	wildcard_constraints:
		size="\d+",
		iter="\d+",
		rep="\d+"
	params:
		output_dir='{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/',
		rep="{rep}"
	shell:
		'mothur "#set.seed(seed={params.rep}); '
		'set.dir(output={params.output_dir}); '
		'cluster(column={input.column}, count={input.count})"'

rule fit:
	input:
		reflist='{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/reference.opti_mcc.list',
		refcolumn='{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/reference.dist',
		refcount='{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/reference.count_table',
		reffasta='{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/reference.fasta',
		fasta='{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/sample.fasta',
		count='{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/sample.count_table',
		column='{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/sample.dist',
	output:
		expand('{{output_dir}}/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{{weight}}_size-{{size}}_i-{{iter}}/r-{{rep}}/method-{{method}}_printref-{{printref}}/sample.opti_mcc.{ext}', ext={'list', 'steps', 'sensspec'})
	wildcard_constraints:
		size="\d+",
		iter="\d+",
		rep="\d+"
	params:
		input_dir="{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/",
		output_dir="{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/method-{method}_printref-{printref}",
		rep="{rep}",
		method="{method}",
		printref='{printref}'
	shell:
		'mothur "set.seed(seed={params.rep}); '
		'set.dir(input={params.input_dir}, output={params.output_dir}); '
		'cluster.fit(reflist={input.reflist}, refcolumn={input.refcolumn}, refcount={input.refcount}, reffasta={input.reffasta}, fasta={input.fasta}, count={input.count}, column={input.column}, printref={params.printref}, method={params.method})'
"""
rule aggregate_sensspec:
	input:
		opticlust=expand('{{output_dir}}/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/{prefix}.opti_mcc.sensspec', weight=weights, size=sizes, iter=iters, rep=reps, prefix=['sample', 'reference']),
		optifit=expand('{{output_dir}}/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/sample.opti_mcc.sensspec', weight=weights, size=sizes, iter=iters, rep=reps, method=methods, printref=printrefs)
	output:
		"{output_dir}/dataset-as-reference/{dataset}/aggregate.sensspec"
	params:
		output_dir="{output_dir}",
		dataset="{dataset}",
		sizes=sizes,
		weights=weights,
		iters=iters,
		reps=reps,
		methods=methods,
		printrefs=printrefs,
		prefixes=['sample','reference']
	shell:
		"for weight in {params.weights}; do "
		"	for size in {params.sizes}; do "
		"		for iter in {params.iters}; do "
		"			for rep in {params.reps}; do "
		"				for prefix in {params.prefixes}; do "
		"					result=$(head -2 ) "
		REF=$(head -2 ${OUTPUTDIR}${REFPI}_${I}_${J}/${PREFIX}reference.opti_mcc.sensspec | tail -1 | sed 's_\(\S*\t\S*\t\)\(.*\)_\t\1\t\2_')
		SAMP=$(head -2 ${OUTPUTDIR}${REFPI}_${I}_${J}/${PREFIX}sample.opti_mcc.sensspec | tail -1 | sed 's_\(\S*\t\S*\t\)\(.*\)_\t\1\t\2_')
		SAMP_O_REF=$(head -2 ${OUTPUTDIR}${REFPI}_${I}_${J}/${PREFIX}sample.open.ref.sensspec | tail -1)
		SAMP_C_REF=$(head -2 ${OUTPUTDIR}${REFPI}_${I}_${J}/${PREFIX}sample.closed.ref.sensspec | tail -1)
		SAMP_O_NOREF=$(head -2 ${OUTPUTDIR}${REFPI}_${I}_${J}/${PREFIX}sample.open.noref.sensspec | tail -1)
		SAMP_C_NOREF=$(head -2 ${OUTPUTDIR}${REFPI}_${I}_${J}/${PREFIX}sample.closed.noref.sensspec | tail -1)

		#REF and SAMP were run with opticlust, which produces sensspec files with 2 less columns than optifit
		#Add two extra tabs at the beginning of their lines so that confusion matrix values line up
		#REF and SAMP also have records that end in a tab, so one less tab at the end
		echo "${REF}${REFP}	$I	$J	REF" >> $FINAL
		echo "${SAMP}${REFP}	$I	$J	SAMP" >> $FINAL
		echo "$SAMP_O_REF	$REFP	$I	$J	SAMP_O_REF" >> $FINAL
		echo "$SAMP_C_REF	$REFP	$I	$J	SAMP_C_REF" >> $FINAL
		echo "$SAMP_O_NOREF	$REFP	$I	$J	SAMP_O_NOREF" >> $FINAL
		echo "$SAMP_C_NOREF	$REFP	$I	$J	SAMP_C_NOREF" >> $FINAL
"""
