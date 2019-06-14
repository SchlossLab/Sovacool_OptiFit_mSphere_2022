""" Benchmarking the OptiFit algorithm """

import math
import os
import re
import shutil

datasets = [dataset_name if not config['subsample_test'] else "{}_{}".format(dataset_name, config['subsample_size']) for dataset_name in config['datasets']]

rule fit_all:
    input:
        expand("results/{output_dir}/{dataset}/figures/aggregate.sensspec.mcc{suffix}.png", output_dir=output_dirs, dataset=datasets, suffix={'', '.full', '.iters'}),
        expand('results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/{prefix}.opti_mcc.sensspec', dataset=datasets, weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, prefix=['sample', 'reference']),
        expand('results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/sample.optifit_mcc.sensspec', dataset=datasets, weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, method=methods, printref=printrefs)

rule split_weighted_subsample:  # TODO: use mothur refweight option instead of this rule
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
        "results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample/sample.dist",
        "results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample/sample.accnos"
    params:
        mothur=mothur_bin,
        output_dir="results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample/",
        iter="{iter}"
    benchmark:
        "benchmarks/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/prep_weighted_subsample.log"
    log:
        "logfiles/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/prep_weighted_subsample.log"
    shell:
        '{params.mothur} "#set.logfile(name={log}); set.seed(seed={params.iter}); set.dir(output={params.output_dir}); get.seqs(accnos={input.accnos}, fasta={input.fasta}); get.seqs(accnos={input.accnos}, count={input.count}); get.dists(column={input.dist}, accnos=current); rename.file(fasta=current, count=current, accnos=current, column=current, prefix=sample)"'

rule prep_reference_from_dataset:
    input:
        fasta=f"{input_dir}/{{dataset}}/{{dataset}}.fasta",
        count=f"{input_dir}/{{dataset}}/{{dataset}}.count_table",
        dist=f"{input_dir}/{{dataset}}/{{dataset}}.dist",
        accnos="results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/sample/sample.accnos"
    output:
        "results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference/reference.accnos",
        "results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference/reference.count_table",
        "results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference/reference.fasta",
        "results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference/reference.dist"
    params:
        mothur=mothur_bin,
        output_dir="results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference/",
        iter="{iter}"
    benchmark:
        "benchmarks/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/prep_reference_from_dataset.log"
    log:
        "logfiles/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/prep_reference_from_dataset.log"
    shell:
        '{params.mothur} "#set.logfile(name={log}); set.seed(seed={params.iter}); set.dir(output={params.output_dir}); remove.seqs(fasta={input.fasta}, count={input.count}, accnos={input.accnos}); list.seqs(fasta=current); get.dists(column={input.dist}, accnos=current); rename.file(fasta=current, count=current, column=current, accnos=current, prefix=reference)"'

rule cluster:
    input:
        count="results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/{sampleref}/{sampleref}.count_table",
        column="results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/{sampleref}/{sampleref}.dist"
    output:
        expand('results/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{{weight}}_reference-fraction-{{reference_fraction}}_i-{{iter}}/r-{{rep}}/{{sampleref}}.opti_mcc.{ext}', ext={'list', 'steps', 'sensspec'})
    params:
        mothur=mothur_bin,
        output_dir='results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/',
        rep="{rep}"
    benchmark:
        "benchmarks/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/{sampleref}.cluster.log"
    log:
        "logfiles/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/{sampleref}.cluster.log"
    shell:
        '{params.mothur} "#set.logfile(name={log}); set.seed(seed={params.rep}); set.dir(output={params.output_dir}); cluster(column={input.column}, count={input.count}, cutoff=0.3)"'

rule fit_to_self:
    input:
        count=f"{input_dir}/{{dataset}}/{{dataset}}.count_table",
        column=f"{input_dir}/{{dataset}}/{{dataset}}.dist",
	accnos=f"{input_dir}/{{dataset}}/{{dataset}}.accnos",
	refaccnos="results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/reference/reference.accnos"
    output:
        temp("results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/{dataset}.dist"),
        expand('results/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{{weight}}_reference-fraction-{{reference_fraction}}_i-{{iter}}/r-{{rep}}/method-{{method}}_printref-{{printref}}/{{dataset}}.optifit_mcc.{ext}', ext={'list', 'sensspec'})
    params:
        mothur=mothur_bin,
        output_dir="results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/",
        rep="{rep}",
        method="{method}",
        printref='{printref}',
	dataset="{dataset}"
    benchmark:
        "benchmarks/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/fit.log"
    log:
        "logfiles/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/fit.log"
    shell:
        '{params.mothur} "#set.logfile(name={log}); set.seed(seed={params.rep}); set.dir(output={params.output_dir}); get.dists(column={input.column}, accnos={input.accnos}); rename.file(column=current, prefix={params.dataset});cluster.fit(column=current, count={input.count}, accnos={input.refaccnos}, printref={params.printref}, method={params.method})"'

rule aggregate_sensspec:
    input:
        expand("results/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/{{dataset}}.dist", weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, method=methods, printref=printrefs),
        opticlust=expand('results/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/{prefix}.opti_mcc.sensspec', weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, prefix=['sample', 'reference']),
        optifit=expand('results/dataset-as-reference/{{dataset}}/{{dataset}}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/{{dataset}}.optifit_mcc.sensspec', weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, method=methods, printref=printrefs)
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
            for weight in params.weights:  # TODO: use regex instead of so many nested for loops
                for reference_fraction in params.reference_fractions:
                    for iter in params.iters:
                        for rep in params.reps:
                            for prefix in params.prefixes:
                                input_filename = f'results/dataset-as-reference/{params.dataset}/{params.dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/{prefix}.opti_mcc.sensspec'
                                with open(input_filename, 'r') as input_file:
                                    for line in input_file:
                                        pass
                                    opticlust_result = re.sub(r'(\S*\t\S*\t)(.*)', r'\t\1\t\2', line).rstrip()
                                    output_file.write(f"{opticlust_result}\t{reference_fraction}\t{iter}\t{rep}\t{prefix}\n")
                            for method in params.methods:
                                for printref in params.printrefs:
                                    input_filename = f"results/dataset-as-reference/{params.dataset}/{params.dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/{params.dataset}.optifit_mcc.sensspec"
                                    with open(input_filename, 'r') as input_file:
                                        for line in input_file:
                                            pass
                                        line = line.strip()
                                        output_file.write(f"{line}\t{reference_fraction}\t{iter}\t{rep}\tmethod-{method}_printref-{printref}\n")

rule plot_sensspec:
    input:
        "results/{output_dir}/{dataset}/aggregate.sensspec"
    output:
        combo_mcc="results/{output_dir}/{dataset}/figures/aggregate.sensspec.mcc.png",
        mcc_full="results/{output_dir}/{dataset}/figures/aggregate.sensspec.mcc.full.png",
        iters="results/{output_dir}/{dataset}/figures/aggregate.sensspec.mcc.iters.png"
    benchmark:
        "benchmarks/{output_dir}/{dataset}/plot_sensspec.log"
    script:
        "plot_sensspec.R"

rule fraction_mapped:
    input:
        mapped=sorted(expand("results/{{dataset}}/{{dataset}}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-closed_printref-f/sample.optifit_mcc.list", weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps)),
        original=sorted(expand("results/{{dataset}}/{{dataset}}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-closed_printref-f/sample.count_table", weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps))
    output:
        "results/{dataset}/{dataset}_fraction_mapped.tsv"
    run:
        if len(input.mapped) != len(input.count_table):
            raise ValueError("Unequal number of optifit_mcc.list and count_table files")
        with open(output[0], 'w') as output_file:
            output_file.write('count_table_filename\tmapped_filename\tfraction_mapped\n')
            for mapped_filename, count_table_filename in zip(input.mapped, input.count_table):
                with open(count_table_filename, 'r') as input_file:
                    line = next(input_file)  # first column of all lines except first line
                    input_samples = set([line.split()[0] for line in input_file])
                with open(mapped_filename, 'r') as mapped_file:
                    line = next(mapped_file)
                    line = next(mapped_file) # third column onward of second line, each seq in each OTU delimited by comma
                    mapped_samples = set(seq for column in line.split()[2:] for seq in column.split(','))
                fraction_mapped = len(input_samples.intersection(mapped_samples)) / len(input_samples)
                output_file.write(f'{count_table_filename}\t{mapped_filename}\t{fraction_mapped}\n')
