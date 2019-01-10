" Download references, download & process data, and run tests to benchmark the OptiFit algorithm"

configfile: 'config/config.yaml'

mothur_bin=config['mothur_bin']
input_dir = config['input_dir']
datasets = [dataset_name if not config['subsample_test'] else "{}_{}".format(dataset_name, config['subsample_size']) for dataset_name in config['datasets']]
weights = config['weights']
iters = range(config['iterations'])
reps = range(config['replicates'])
methods = config['methods']
printrefs = config['printrefs']
reference_fractions = [i/100 for i in range(config['reference_fractions']['start'], config['reference_fractions']['stop'], config['reference_fractions']['step'])]
output_dirs = [option for option in config['workflows'] if config['workflows'][option]]

wildcard_constraints:
    dataset="\w+",
    iter="\d+",
    rep="\d+",
    sampleref="sample|reference"


include: 'code/data_processing/get-references.smk'
# TODO: write & include snakemake workflows to replace {dataset}.batch and {dataset}.R files
include: 'code/data_processing/testset-subsample.smk'
include: 'code/analysis/optifit-dataset-as-ref.smk'
#include: 'code/analysis/optifit-silva-ref.smk'


rule all:
        input:
            expand("results/{output_dir}/{dataset}/figures/aggregate.sensspec.mcc{suffix}.png", output_dir=output_dirs, dataset=datasets, suffix={'', '.full', '.iters'})

rule calc_seq_dists:
    input:
        f'{input_dir}/{{dataset}}/{{dataset}}.fasta'
    output:
        f'{input_dir}/{{dataset}}/{{dataset}}.dist'
    params:
        mothur=mothur_bin,
        output_dir=f'{input_dir}/{{dataset}}/'
    benchmark:
        f'benchmarks/{input_dir}/{{dataset}}/calc_seq_dists.log'
    log:
        f"logfiles/{input_dir}/{{dataset}}/calc_seq_dists.log"
    shell:
        '{params.mothur} "#set.logfile(name={log}); set.dir(output={params.output_dir}); dist.seqs(fasta={input[0]}, cutoff=0.03)"'
"""
rule fraction_mapped:
    input:
        mapped=sorted(expand("results/{{output_dir}}/{{dataset}}/{{dataset}}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-closed_printref-f/sample.optifit_mcc.list", weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps)),
        original=sorted(expand("results/{{output_dir}}/{{dataset}}/{{dataset}}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-closed_printref-f/sample.count_table", weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps))
    output:
        "results/{output_dir}/{dataset}/{dataset}_fraction_mapped.tsv"
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
"""
