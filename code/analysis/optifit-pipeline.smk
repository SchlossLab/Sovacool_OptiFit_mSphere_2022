""" Benchmarking the OptiFit algorithm """

import math
import os

configfile: 'config_test.yaml'

input_dir = config['input_dir']
datasets = [dataset_name if not config['subsample_test'] else "{}_{}".format(dataset_name, config['subsample_size']) for dataset_name in config['datasets']]
weights = config['weights']
iters = range(config['iterations'])
reps = range(config['replicates'])
methods = {'open', 'closed'}
printrefs = {'t', 'f'}
reference_fractions = [i/100 for i in range(50,70,10)]

wildcard_constraints:
	dataset="\w+",
	iter="\d+",
	rep="\d+",
	sampleref="\w+"

rule all:
	input:
		expand('results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/{sampleref}.opti_mcc.{ext}', dataset=datasets, weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, sampleref=['sample', 'reference'], ext={'list', 'steps', 'sensspec'}),
		expand('results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/sample.optifit_mcc.{ext}', dataset=datasets, weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, method=methods, printref=printrefs, ext={'list', 'steps', 'sensspec'}),
		expand("results/dataset-as-reference/{dataset}/figures/aggregate.sensspec.mcc{suffix}.png", dataset=datasets, suffix={'', '.full', '.iters'})

rule get_dists:
	input:
		f'{input_dir}/{{dataset}}/{{dataset}}.fasta'
	output:
		temp(f'{input_dir}/{{dataset}}/{{dataset}}.dist')
	params:
		output_dir=f'{input_dir}/{{dataset}}/'
	benchmark:
		f'benchmarks/{input_dir}/{{dataset}}/get_dists.log'
	log:
		f"logfiles/{input_dir}/{{dataset}}/get_dists.log"
	shell:
		'mothur "#set.logfile(name={log}); set.dir(output={params.output_dir}); dist.seqs(fasta={input[0]}, cutoff=0.03)"'

rule split_weighted_subsample:
	input:
		count=f"{input_dir}/{{dataset}}/{{dataset}}.count_table",
		dist=f"{input_dir}/{{dataset}}/{{dataset}}.dist"
	output:
		"results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample.accnos"
	params:
		dataset="{dataset}",
		reference_fraction="{reference_fraction}",
		weight="{weight}"
	benchmark:
		'benchmarks/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/split_weighted_subsample.log'
	script:
		"weighted_subsample.R"

rule prep_weighted_subsample:
	input:
		fasta=f"{input_dir}/{{dataset}}/{{dataset}}.fasta",
		count=f"{input_dir}/{{dataset}}/{{dataset}}.count_table",
		dist=f"{input_dir}/{{dataset}}/{{dataset}}.dist",
		accnos="results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample.accnos"
	output:
		"results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample/sample.fasta",
		"results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample/sample.count_table",
		temp("results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample/sample.dist")
	params:
		output_dir="results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample/",
		iter="{iter}"
	benchmark:
		"benchmarks/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/prep_weighted_subsample.log"
	log:
		"logfiles/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/prep_weighted_subsample.log"
	shell:
		'mothur "#set.logfile(name={log}); set.seed(seed={params.iter}); set.dir(output={params.output_dir}); get.seqs(accnos={input.accnos}, fasta={input.fasta}); get.seqs(accnos={input.accnos}, count={input.count}); get.dists(column={input.dist}, accnos=current); rename.file(fasta=current, count=current, accnos=current, column=current, prefix=sample)"'

rule prep_reference_from_dataset:
	input:
		fasta=f"{input_dir}/{{dataset}}/{{dataset}}.fasta",
		count=f"{input_dir}/{{dataset}}/{{dataset}}.count_table",
		dist=f"{input_dir}/{{dataset}}/{{dataset}}.dist",
		accnos="results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample.accnos"
	output:
		"results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference/reference.accnos",
		"results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference/reference.count_table",
		"results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference/reference.fasta",
		temp("results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference/reference.dist")
	params:
		output_dir="results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference/",
		iter="{iter}"
	benchmark:
		"benchmarks/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/prep_reference_from_dataset.log"
	log:
		"logfiles/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/prep_reference_from_dataset.log"
	shell:
		'mothur "#set.logfile(name={log}); set.seed(seed={params.iter}); set.dir(output={params.output_dir}); remove.seqs(fasta={input.fasta}, count={input.count}, accnos={input.accnos}); list.seqs(fasta=current); get.dists(column={input.dist}, accnos=current); rename.file(fasta=current, count=current, column=current, accnos=current, prefix=reference)"'

rule cluster:
	input:
		count="results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/{sampleref}/{sampleref}.count_table",
		column="results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/{sampleref}/{sampleref}.dist"
	output:
		expand('results/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{{weight}}_reference-fraction-{{reference_fraction}}_i-{{iter}}/r-{{rep}}/{{sampleref}}.opti_mcc.{ext}', ext={'list', 'steps', 'sensspec'})
	params:
		output_dir='results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/',
		rep="{rep}"
	benchmark:
		"benchmarks/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/{sampleref}.cluster.log"
	log:
		"logfiles/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/{sampleref}.cluster.log"
	shell:
		'mothur "#set.logfile(name={log}); set.seed(seed={params.rep}); set.dir(output={params.output_dir}); cluster(column={input.column}, count={input.count}, cutoff=0.3)"'

rule fit:
	input:
		reflist='results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/reference.opti_mcc.list',
		refcolumn='results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference/reference.dist',
		refcount='results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference/reference.count_table',
		reffasta='results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference/reference.fasta',
		fasta='results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample/sample.fasta',
		count='results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample/sample.count_table',
		column='results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample/sample.dist'
	output:
		expand('results/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{{weight}}_reference-fraction-{{reference_fraction}}_i-{{iter}}/r-{{rep}}/method-{{method}}_printref-{{printref}}/sample.optifit_mcc.{ext}', ext={'list', 'steps', 'sensspec'}),
		temp(expand('results/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{{weight}}_reference-fraction-{{reference_fraction}}_i-{{iter}}/r-{{rep}}/method-{{method}}_printref-{{printref}}/sample.{suffix}.dist', suffix={'pick', 'fit'}))
	params:
		#input_dir="results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/",
		output_dir="results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/",
		rep="{rep}",
		method="{method}",
		printref='{printref}'
	benchmark:
		"benchmarks/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/fit.log"
	log:
		"logfiles/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/fit.log"
	shell:
		'mothur "#set.logfile(name={log}); set.seed(seed={params.rep}); set.dir(output={params.output_dir}); cluster.fit(reflist={input.reflist}, refcolumn={input.refcolumn}, refcount={input.refcount}, reffasta={input.reffasta}, fasta={input.fasta}, count={input.count}, column={input.column}, printref={params.printref}, method={params.method})"'

rule aggregate_sensspec:
	input:
		opticlust=expand('results/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/{prefix}.opti_mcc.sensspec', weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, prefix=['sample', 'reference']),
		optifit=expand('results/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/sample.optifit_mcc.sensspec', weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, method=methods, printref=printrefs)
	output:
		"results/dataset-as-reference/{dataset}/aggregate.sensspec"
	params:
		dataset="{dataset}",
		reference_fractions=reference_fractions,
		weights=weights,
		iters=iters,
		reps=reps,
		methods=methods,
		printrefs=printrefs,
		prefixes=['sample','reference']
	benchmark:
		"benchmarks/dataset-as-reference/{dataset}/aggregate_sensspec.log"
	shell:
		"if [ -e {output[0]} ]; then rm {output[0]} ; fi ; "
		"echo 'iter	label	cutoff	numotus	tp	tn	fp	fn	sensitivity	specificity	ppv	npv	fdr	accuracy	mcc	f1score	reference_fraction	iter	rep	type' >> {output[0]} "
		"reference_fractions=({params.reference_fractions}); weights=({params.weights}); iters=({params.iters}); reps=({params.reps}); methods=({params.methods}); printrefs=({params.printrefs}); prefixes=({params.prefixes}); "
		"for weight in ${{weights[@]}}; do "
		"	for reference_fraction in ${{reference_fractions[@]}}; do "
		"		for iter in ${{iters[@]}}; do "
		"			for rep in ${{reps[@]}}; do "
		"				for prefix in ${{prefixes[@]}}; do "
		"					opticlust=$(head -2 results/dataset-as-reference/{params.dataset}/{params.dataset}_weight-${{weight}}_reference-fraction-${{reference_fraction}}_i-${{iter}}/r-${{rep}}/${{prefix}}.opti_mcc.sensspec | tail -1 | sed 's_\(\S*\t\S*\t\)\(.*\)_\t\1\t\2_' ); "
		"					echo '${{opticlust}}${{reference_fraction}} ${{iter}} ${{rep}} ${{prefix}}' >> {output[0]}; "
		"				done; "
		"				for method in ${{methods[@]}}; do "
		"					for printref in ${{params.printrefs[@]}}; do "
		"						optifit=$(head -2 results/dataset-as-reference/{params.dataset}/{params.dataset}_weight-${{weight}}_reference-fraction-${{reference_fraction}}_i-${{iter}}/r-${{rep}}/method-${{method}}_printref-${{printref}}/sample.opti_mcc.sensspec | tail -1); "
		"						echo '${{optifit}} ${{reference_fraction}} ${{iter}} ${{rep}} method-${{method}}_printref-${{printref}}' >> {output[0]}; "
		"					done; "
		"				done; "
		"			done; "
		"		done; "
		"	done: "
		"done; "

rule plot_sensspec:
	input:
		"results/dataset-as-reference/{dataset}/aggregate.sensspec"
	output:
		combo_mcc="results/dataset-as-reference/{dataset}/figures/aggregate.sensspec.mcc.png",
		mcc_full="results/dataset-as-reference/{dataset}/figures/aggregate.sensspec.mcc.full.png",
		iters="results/dataset-as-reference/{dataset}/figures/aggregate.sensspec.mcc.iters.png"
	benchmark:
		"benchmarks/dataset-as-reference/{dataset}/plot_sensspec.log"
	script:
		"plot_sensspec.R"
