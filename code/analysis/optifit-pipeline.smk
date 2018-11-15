""" Benchmarking the OptiFit algorithm """

import math
import os

configfile: 'config_test.yaml'

input_dir = config['input_dir']
output_dir = config['output_dir']
datasets = [dataset_name if not config['subsample_test'] else "{}_{}".format(dataset_name, config['subsample_size']) for dataset_name in config['datasets']]
weights = config['weights']
iters = range(config['iterations'])
reps = range(config['replicates'])
methods = {'open', 'closed'}
printrefs = {'t', 'f'}
reference_fractions = [i/100 for i in range(50,70,10)]

rule all:
	input:
		expand('{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/{sampleref}.opti_mcc.{ext}', output_dir=output_dir, dataset=datasets, weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, sampleref=['sample', 'reference'], ext={'list', 'steps', 'sensspec'}),
		expand('{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/sample.opti_mcc.{ext}', output_dir=output_dir, dataset=datasets, weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, method=methods, printref=printrefs, ext={'list', 'steps', 'sensspec'}),
		expand("{output_dir}/dataset-as-reference/{dataset}/figures/aggregate.sensspec.mcc{suffix}.png", output_dir=output_dir, dataset=datasets, suffix={'', '.full', '.iters'})

rule get_dists:
	input:
		'{input_dir}/{dataset}/{dataset}.fasta'
	output:
		temp('{input_dir}/{dataset}/{dataset}.dist')
	params:
		output_dir='{input_dir}/{dataset}/'
	shell:
		'mothur "#set.dir(output={params.output_dir}); dist.seqs(fasta={input[0]}, cutoff=0.03)"'

rule split_weighted_subsample:
	input:
		count=f"{input_dir}/{{dataset}}/{{dataset}}.count_table",
		dist=f"{input_dir}/{{dataset}}/{{dataset}}.dist"
	output:
		"{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample.accnos"
	params:
		dataset="{dataset}",
		reference_fraction="{reference_fraction}",
		weight="{weight}"
	wildcard_constraints:
		iter="\d+"
	script:
		"weighted_subsample.R"

rule prep_weighted_subsample:
	input:
		fasta=f"{input_dir}/{{dataset}}/{{dataset}}.fasta",
		count=f"{input_dir}/{{dataset}}/{{dataset}}.count_table",
		dist=f"{input_dir}/{{dataset}}/{{dataset}}.dist",
		accnos="{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample.accnos"
	output:
		"{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample.fasta",
		"{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample.count_table",
		temp("{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample.dist")
	params:
		output_dir="{output_dir}/{dataset}/dataset-as-reference/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/",
		iter="{iter}"
	wildcard_constraints:
		iter="\d+"
	shell:
		'mothur "#set.seed(seed={params.iter}); set.dir(output={params.output_dir}); get.seqs(accnos={input.accnos}, fasta={input.fasta}, count={input.count}); get.dists(column={input.dist}, accnos=current); rename.file(fasta=current, count=current, accnos = current, column=current, prefix=sample)"'

rule prep_reference_from_dataset:
	input:
		fasta="{input_dir}/{{dataset}}/{{dataset}}.fasta".format(input_dir=input_dir),
		count="{input_dir}/{{dataset}}/{{dataset}}.count_table".format(input_dir=input_dir),
		dist="{input_dir}/{{dataset}}/{{dataset}}.dist".format(input_dir=input_dir),
		accnos="{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample.accnos"
	output:
		"{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference.accnos",
		"{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference.count_table",
		"{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference.fasta",
		temp("{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference.dist")
	wildcard_constraints:
		iter="\d+"
	params:
		output_dir="{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/",
		iter="{iter}"
	shell:
		'mothur "#set.seed(seed={params.iter}); set.dir(output={params.output_dir}); remove.seqs(fasta={input.fasta}, count={input.count}, accnos={input.accnos}); list.seqs(fasta=current); get.dists(column={input.dist}, accnos=current); rename.file(fasta=current, count=current, column=current, accnos=current, prefix=reference)"'

rule cluster:
	input:
		count="{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/{sampleref}.count_table",
		column="{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/{sampleref}.dist"
	output:
		expand('{{output_dir}}/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{{weight}}_reference-fraction-{{reference_fraction}}_i-{{iter}}/r-{{rep}}/{{sampleref}}.opti_mcc.{ext}', ext={'list', 'steps', 'sensspec'})
	wildcard_constraints:
		iter="\d+",
		rep="\d+",
		sampleref="\w+"
	params:
		output_dir='{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/',
		rep="{rep}"
	shell:
		'mothur "#set.seed(seed={params.rep}); set.dir(output={params.output_dir}); cluster(column={input.column}, count={input.count})"'

rule fit:
	input:
		reflist='{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/reference.opti_mcc.list',
		refcolumn='{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference.dist',
		refcount='{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference.count_table',
		reffasta='{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference.fasta',
		fasta='{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample.fasta',
		count='{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample.count_table',
		column='{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample.dist'
	output:
		expand('{{output_dir}}/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{{weight}}_reference-fraction-{{reference_fraction}}_i-{{iter}}/r-{{rep}}/method-{{method}}_printref-{{printref}}/sample.opti_mcc.{ext}', ext={'list', 'steps', 'sensspec'}),
		temp(expand('{{output_dir}}/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{{weight}}_reference-fraction-{{reference_fraction}}_i-{{iter}}/r-{{rep}}/method-{{method}}_printref-{{printref}}/sample.{suffix}.dist', suffix={'pick', 'fit'}))
	params:
		input_dir="{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/",
		output_dir="{output_dir}/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}",
		rep="{rep}",
		method="{method}",
		printref='{printref}'
	wildcard_constraints:
		iter="\d+",
		rep="\d+"
	shell:
		'mothur "set.seed(seed={params.rep}); set.dir(input={params.input_dir}, output={params.output_dir}); cluster.fit(reflist={input.reflist}, refcolumn={input.refcolumn}, refcount={input.refcount}, reffasta={input.reffasta}, fasta={input.fasta}, count={input.count}, column={input.column}, printref={params.printref}, method={params.method})'

rule aggregate_sensspec:
	input:
		opticlust=expand('{{output_dir}}/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/{prefix}.opti_mcc.sensspec', weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, prefix=['sample', 'reference']),
		optifit=expand('{{output_dir}}/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/sample.opti_mcc.sensspec', weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, method=methods, printref=printrefs)
	output:
		"{output_dir}/dataset-as-reference/{dataset}/aggregate.sensspec"
	params:
		output_dir="{output_dir}",
		dataset="{dataset}",
		reference_fractions=reference_fractions,
		weights=weights,
		iters=iters,
		reps=reps,
		methods=methods,
		printrefs=printrefs,
		prefixes=['sample','reference']
	wildcard_constraints:
		iter="\d+",
		rep="\d+"
	shell:
		"if [ -e {output[0]} ]; then rm {output[0]} ; fi ; "
		"echo 'iter	label	cutoff	numotus	tp	tn	fp	fn	sensitivity	specificity	ppv	npv	fdr	accuracy	mcc	f1score	reference_fraction	iter	rep	type' >> {output[0]} "
		"reference_fractions=({params.reference_fractions}); weights=({params.weights}); iters=({params.iters}); reps=({params.reps}); methods=({params.methods}); printrefs=({params.printrefs}); prefixes=({params.prefixes}); "
		"for weight in ${{weights[@]}}; do "
		"	for reference_fraction in ${{reference_fractions[@]}}; do "
		"		for iter in ${{iters[@]}}; do "
		"			for rep in ${{reps[@]}}; do "
		"				for prefix in ${{prefixes[@]}}; do "
		"					opticlust=$(head -2 {params.output_dir}/dataset-as-reference/{params.dataset}/{params.dataset}_weight-${{weight}}_reference-fraction-${{reference_fraction}}_i-${{iter}}/r-${{rep}}/${{prefix}}.opti_mcc.sensspec | tail -1 | sed 's_\(\S*\t\S*\t\)\(.*\)_\t\1\t\2_' ); "
		"					echo '${{opticlust}}${{reference_fraction}} ${{iter}} ${{rep}} ${{prefix}}' >> {output[0]}; "
		"				done; "
		"				for method in ${{methods[@]}}; do "
		"					for printref in ${{params.printrefs[@]}}; do "
		"						optifit=$(head -2 {params.output_dir}/dataset-as-reference/{params.dataset}/{params.dataset}_weight-${{weight}}_reference-fraction-${{reference_fraction}}_i-${{iter}}/r-${{rep}}/method-${{method}}_printref-${{printref}}/sample.opti_mcc.sensspec | tail -1); "
		"						echo '${{optifit}} ${{reference_fraction}} ${{iter}} ${{rep}} method-${{method}}_printref-${{printref}}' >> {output[0]}; "
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
