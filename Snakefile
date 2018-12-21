" Download references, download & process data, and run tests to benchmark the OptiFit algorithm"

configfile: 'config.yaml'

mothur_bin=config['mothur_bin']
input_dir = config['input_dir']
datasets = [dataset_name if not config['subsample_test'] else "{}_{}".format(dataset_name, config['subsample_size']) for dataset_name in config['datasets']]
weights = config['weights']
iters = range(config['iterations'])
reps = range(config['replicates'])
methods = {'open', 'closed'}
printrefs = {'t', 'f'}
reference_fractions = [i/100 for i in range(50,60,10)]
output_dirs = [option for option in config['workflows'] if config['workflows'][option]]


wildcard_constraints:
    dataset="\w+",
    iter="\d+",
    rep="\d+",
    sampleref="sample|reference",
    reference="silva|greengenes"


include: 'code/data_processing/get-references.smk'
# TODO: write & include snakemake workflows to replace {dataset}.batch and {dataset}.R files
include: 'code/data_processing/testset-subsample.smk'
include: 'code/analysis/optifit-dataset-as-ref.smk'
include: 'code/analysis/optifit-ref-db.smk'

rule all:
        input:
            expand("results/{output_dir}/{dataset}/figures/aggregate.sensspec.mcc{suffix}.png", output_dir=output_dirs, dataset=datasets, suffix={'', '.full', '.iters'})

rule calc_seq_dists:
    input:
        '{input_dir}/{sample}/{sample}.fasta'
    output:
        '{input_dir}/{sample}/{sample}.dist'
    params:
        mothur=mothur_bin,
        output_dir='{input_dir}/{sample}/'
    benchmark:
        'benchmarks/{input_dir}/{sample}/calc_seq_dists.log'
    log:
        "logfiles/{input_dir}/{sample}/calc_seq_dists.log"
    shell:
        '{params.mothur} "#set.logfile(name={log}); set.dir(output={params.output_dir}); dist.seqs(fasta={input[0]}, cutoff=0.03)"'

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
