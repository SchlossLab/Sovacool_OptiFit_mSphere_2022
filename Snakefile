" Download references, download & process data, and run tests to benchmark the OptiFit algorithm"

configfile: 'config_test.yaml'
include: 'code/data_processing/get-references.smk'
# TODO: write & include snakemake workflows to replace {dataset}.batch and {dataset}.R files
include: 'code/data_processing/testset-subsample.smk'
include: 'code/analysis/optifit-pipeline.smk'
