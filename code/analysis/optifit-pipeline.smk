""" Benchmarking the OptiFit algorithm """

import math
import os
import re

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
		f'{input_dir}/{{dataset}}/{{dataset}}.dist'
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
		"results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample/sample.dist"
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
		"results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference/reference.dist"
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
		expand('results/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{{weight}}_reference-fraction-{{reference_fraction}}_i-{{iter}}/r-{{rep}}/method-{{method}}_printref-{{printref}}/sample.optifit_mcc.{ext}', ext={'list', 'steps', 'sensspec'})
	params:
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
	run:
		header_str = 'iter\tlabel\tcutoff\tnumotus\ttp\ttn\tfp\tfn\tsensitivity\tspecificity\tppv\tnpv\tfdr\taccuracy\tmcc\tf1score\treference_fraction\titer\trep\ttype\n'
		with open(output[0], 'w') as output_file:
			output_file.write(header_str)
			for weight in params.weights:
				for reference_fraction in params.reference_fractions:
					for iter in params.iters:
						for rep in params.reps:
							for prefix in params.prefixes:
								input_filename = f'results/dataset-as-reference/{params.dataset}/{params.dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/{prefix}.opti_mcc.sensspec'
								with open(input_filename, 'r') as input_file:
									for line in input_file:
										pass
									opticlust_result = re.sub("\(\S*\t\S*\t\)\(.*\)", "\t\1\t\2", line.strip())
									output_file.write(f"{opticlust_result}\t{reference_fraction}\t{iter}\t{rep}\t{prefix}\n")
							for method in params.methods:
								for printref in params.printrefs:
									input_filename = f"results/dataset-as-reference/{params.dataset}/{params.dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/sample.optifit_mcc.sensspec"
									with open(input_filename, 'r') as input_file:
										for line in input_file:
											pass
										line = line.strip()
										output_file.write(f"{line}\t{reference_fraction}\t{iter}\t{rep}\tmethod-{method}_printref-{printref}\n")

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
