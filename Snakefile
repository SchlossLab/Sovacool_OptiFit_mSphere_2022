" Download references, download & process results, and run tests to benchmark the OptiFit algorithm "

configfile: 'config/config.yaml'

datasets = [sample_name if not config['subsample_test'] else "{}_{}".format(sample_name, config['subsample_size']) for sample_name in config['datasets']]
weights = set(config['weights'])
methods = set(config['methods'])
printrefs = set(config['printrefs'])
start = int(config['reference_fractions']['start'])
stop = int(config['reference_fractions']['stop'])
step = int(config['reference_fractions']['step'])
reference_fractions = [i/100 for i in range(start, stop, step)]
seeds = range(int(config['seeds']))

subworkflow prep_db:
    workdir:
        "subworkflows/0_prep_db"

subworkflow prep_samples:
    workdir:
        "subworkflows/1_prep_samples"

rule targets:
    input:
        [prep_db(f"results/{ref}/{ref}.{region}.seed_{seed}.opti_mcc.0.03.cons.taxonomy"), for ref in ("silva", "gg", "rdp") for region in ("bact_v4", "bact_full") for seed in seeds],
        [prep_samples(f"results/{dataset}/{dataset}.seed_{seed}.opti_mcc.sensspec") for dataset in datasets for seed in seeds]
