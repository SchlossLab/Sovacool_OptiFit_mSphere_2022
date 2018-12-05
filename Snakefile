" Download references, download & process data, and run tests to benchmark the OptiFit algorithm"

configfile: 'config.yaml'

mothur_bin=config['mothur_bin']
input_dir = config['input_dir']
datasets = [dataset_name if not config['subsample_test'] else "{}_{}".format(dataset_name, config['subsample_size']) for dataset_name in config['datasets']]
weights = config['weights']
iters = range(config['iterations'])
reps = range(config['replicates'])
methods = config['methods']
printrefs = config['printrefs']
reference_fractions = [i/100 for i in range(config['reference_fractions']['start'], config['reference_fractions']['stop'], config['reference_fractions']['step'])]

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

output_dirs = [option for option in ['dataset-as-reference', 'silva-as-reference'] if config[option]]

rule all:
        input:
            expand('results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/{prefix}.opti_mcc.sensspec', dataset=datasets, weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, prefix=['sample', 'reference']),
            expand('results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/sample.optifit_mcc.sensspec', dataset=datasets, weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, method=methods, printref=printrefs),
            expand("results/dataset-as-reference/{dataset}/figures/aggregate.sensspec.mcc{suffix}.png", dataset=datasets, suffix={'', '.full', '.iters'})

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
