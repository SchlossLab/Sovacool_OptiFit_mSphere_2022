" Download references, download & process data, and run tests to benchmark the OptiFit algorithm "

configfile: 'config/config.yaml'

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

subworkflow prep_db:
    workdir:
        "subworkflows/0_prep_db"

rule targets:
    input:
        prep_db("data/silva/silva.bact_v4.filter.unique.precluster.opti_mcc.0.03.cons.taxonomy"),
        prep_db("data/silva/silva.bact_full.filter.unique.precluster.opti_mcc.0.03.cons.taxonomy"),
        prep_db("data/gg/gg.bact_v4.filter.unique.precluster.opti_mcc.0.03.cons.taxonomy"),
        prep_db("data/gg/gg.bact_full.filter.unique.precluster.opti_mcc.0.03.cons.taxonomy"),
        prep_db("data/rdp/rdp.bact_v4.filter.unique.precluster.opti_mcc.0.03.cons.taxonomy"),
        prep_db("data/rdp/rdp.bact_full.filter.unique.precluster.opti_mcc.0.03.cons.taxonomy")
