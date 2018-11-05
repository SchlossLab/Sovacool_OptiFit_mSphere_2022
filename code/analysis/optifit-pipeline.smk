" Benchmarking the OptiFit algorithm "

import os

configfile: 'config_soil.yaml'
dataset = config['dataset']
input_dir = os.path.join(config['input_dir'], dataset)
output_dir = os.path.join(config['output_dir'], dataset)

# TODO: include snakemake workflows to download data (in meantime, call scripts in code/data_processing)
