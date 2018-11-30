" Download references, download & process data, and run tests to benchmark the OptiFit algorithm"

configfile: 'config.yaml'

input_dir = config['input_dir']
datasets = [dataset_name if not config['subsample_test'] else "{}_{}".format(dataset_name, config['subsample_size']) for dataset_name in config['datasets']]
weights = config['weights']
iters = range(config['iterations'])
reps = range(config['replicates'])
methods = {'open', 'closed'}
printrefs = {'t', 'f'}
reference_fractions = [i/100 for i in range(50,60,10)]

wildcard_constraints:
        dataset="\w+",
        iter="\d+",
        rep="\d+",
        sampleref="sample|reference"


include: 'code/data_processing/get-references.smk'
# TODO: write & include snakemake workflows to replace {dataset}.batch and {dataset}.R files
include: 'code/data_processing/testset-subsample.smk'
include: 'code/analysis/optifit-pipeline.smk'

rule all:
        input:
            expand("results/dataset-as-reference/{dataset}/figures/aggregate.sensspec.mcc{suffix}.png", dataset=datasets, suffix={'', '.full', '.iters'})
