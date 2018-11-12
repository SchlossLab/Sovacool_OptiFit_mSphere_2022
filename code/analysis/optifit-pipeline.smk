""" Benchmarking the OptiFit algorithm """

import math
import os

configfile: 'config_test.yaml'

input_dir = config['input_dir']
output_dir = config['output_dir']
datasets = [dataset_name if config['subsample_test'] else "{}_{}".format(dataset_name, config['subsample_size']) for dataset_name in config['datasets']]
weights = config['weights']
iters = range(config['iterations'])
reps = range(config['replicates'])
methods = {'open', 'closed'}
printrefs = {'t', 'f'}
sizes = [math.floor(i/100) for i in range(50,61,10)]

rule all:
	input:
		expand('{input_dir}/{dataset}/{dataset}.{ext}', dataset=datasets, input_dir=input_dir, ext={'fasta', 'dist', 'count_table'}),
		expand("{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.{ext}", output_dir=output_dir, dataset=datasets, weight=weights, size=sizes, iter=iters, ext={'fasta','count_table', 'dist'}),
		expand('{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/sample.opti_mcc.{ext}', output_dir=output_dir, dataset=datasets, weight=weights, size=sizes, iter=iters, rep=reps, method=methods, printref=printrefs, ext={'list', 'steps', 'sensspec'}),
		expand("{output_dir}/dataset-as-reference/{dataset}/figures/aggregate.sensspec.mcc{suffix}.png", output_dir=output_dir, dataset=datasets, suffix={'', '.full', '.iters'})

rule get_dists:
	input:
		'{input_dir}/{dataset}/{dataset}.fasta'
	output:
		'{input_dir}/{dataset}/{dataset}.dist'
	params:
		output_dir='{input_dir}/{dataset}/'
	shell:
		'mothur "#set.dir(output={params.output_dir}); '
		'dist.seqs(fasta={input[0]}, cutoff=0.03);"'

rule split_weighted_subsample:
	input:
		count="{input_dir}/{{dataset}}/{{dataset}}.count_table".format(input_dir=input_dir),
		dist="{input_dir}/{{dataset}}/{{dataset}}.dist".format(input_dir=input_dir)
	output:
		"{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.accnos"
	params:
		dataset="{dataset}",
		size="{size}",
		weight="{weight}"
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
		"if [ -e {output[0]} ]; then rm {output[0]} ; fi ; "
		"echo 'iter	label	cutoff	numotus	tp	tn	fp	fn	sensitivity	specificity	ppv	npv	fdr	accuracy	mcc	f1score	refp	refpi	iter	type' >> {output[0]} "
		"sizes=({params.sizes}); weights=({params.weights}); iters=({params.iters}); reps=({params.reps}); methods=({params.methods}); printrefs=({params.printrefs}); prefixes=({params.prefixes}); "
		"for weight in ${{weights[@]}}; do "
		"	for size in ${{sizes[@]}}; do "
		"		for iter in ${{iters[@]}}; do "
		"			for rep in ${{reps[2]}}; do "
		"				for prefix in ${{prefixes[@]}}; do "
		"					opticlust=$(head -2 {params.output_dir}/dataset-as-reference/{params.dataset}/{params.dataset}_weight-${{weight}}_size-${{size}}_i-${{iter}}/r-${{rep}}/${{prefix}}.opti_mcc.sensspec | tail -1 | sed 's_\(\S*\t\S*\t\)\(.*\)_\t\1\t\2_' ); "
		"					echo '${{opticlust}}${{size}} ${{iter}} ${{rep}} ${{prefix}}' >> {output[0]}; "
		"				done; "
		"				for method in ${{methods[@]}}; do "
		"					for printref in ${{params.printrefs[@]}}; do "
		"						optifit=$(head -2 {params.output_dir}/dataset-as-reference/{params.dataset}/{params.dataset}_weight-${{weight}}_size-${{size}}_i-${{iter}}/r-${{rep}}/method-${{method}}_printref-${{printref}}/sample.opti_mcc.sensspec | tail -1); "
		"						echo '${{optifit}} ${{size}} ${{iter}} ${{rep}} method-${{method}}_printref-${{printref}}' >> {output[0]}; "
		"					done; "
		"				done; "
		"			done; "
		"		done; "
		"	done: "
		"done; "

rule plot_sensspec:
	input:
		"{output_dir}/dataset-as-reference/{dataset}/aggregate.sensspec"
	output:
		combo_mcc="{output_dir}/dataset-as-reference/{dataset}/figures/aggregate.sensspec.mcc.png",
		mcc_full="{output_dir}/dataset-as-reference/{dataset}/figures/aggregate.sensspec.mcc.full.png",
		iters="{output_dir}/dataset-as-reference/{dataset}/figures/aggregate.sensspec.mcc.iters.png"
	script:
		"plot_sensspec.R"
