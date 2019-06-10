" Download references, download & process data, and run tests to benchmark the OptiFit algorithm "

configfile: 'config/config.yaml'

mothur_bin = config['mothur_bin']
input_dir = config['input_dir']
samples = [sample_name if not config['subsample_test'] else "{}_{}".format(sample_name, config['subsample_size']) for sample_name in config['samples']]
weights = set(config['weights'])
methods = set(config['methods'])
printrefs = set(config['printrefs'])
start = int(config['reference_fractions']['start'])
stop = int(config['reference_fractions']['stop'])
step = int(config['reference_fractions']['step'])
reference_fractions = [i/100 for i in range(start, stop, step)]
iters = range(config['iterations'])
reps = range(config['replicates'])

wildcard_constraints:
    sample="\w+",
    iter="\d+",
    rep="\d+",
    sampleref="sample|reference",
    reference="silva|greengenes"
