" Download references, download & process data, and run tests to benchmark the OptiFit algorithm"

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
        sampleref="sample|reference"


include: 'code/data_processing/get-references.smk'
# TODO: write & include snakemake workflows to replace {dataset}.batch and {dataset}.R files
include: 'code/data_processing/testset-subsample.smk'
include: 'code/analysis/optifit-pipeline.smk'

rule all:
        input:
                expand('results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/{sampleref}.opti_mcc.{ext}', dataset=datasets, weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, sampleref=['sample', 'reference'], ext={'list', 'steps', 'sensspec'}),
                expand('results/dataset-as-reference/{dataset}/{dataset}_weight-{weight}_reference-fraction-{reference_fraction}_i-{iter}/r-{rep}/method-{method}_printref-{printref}/sample.optifit_mcc.{ext}', dataset=datasets, weight=weights, reference_fraction=reference_fractions, iter=iters, rep=reps, method=methods, printref=printrefs, ext={'list', 'steps', 'sensspec'}),
                expand("results/dataset-as-reference/{dataset}/figures/aggregate.sensspec.mcc{suffix}.png", dataset=datasets, suffix={'', '.full', '.iters'})


