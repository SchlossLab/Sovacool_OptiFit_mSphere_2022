""" Benchmarking the OptiFit algorithm """

import math
import os

configfile: 'config_soil.yaml'

input_dir = os.path.join(config['input_dir'], config['dataset'])
output_dir = os.path.join(config['output_dir'], config['dataset'])
dataset = config['dataset'] if not config['subsample_test'] else '.'.join([str(config['subsample_size']), config['dataset']])

def count_input_seqs(infilename):
	with open(infilename, 'r') as infile:
		num_seqs = len([line for line in infile if line[0] == '>'])
	return num_seqs
num_seqs = count_input_seqs(os.path.join(input_dir, dataset + '.fasta'))
sizes = [math.floor(num_seqs * i/100) for i in range(5,100,5)]

rule all:
	input:
		expand("{input_dir}/subsamples/size={size}_weight={weight}_i={iter}_r={rep}/sample.accnos", input_dir=input_dir, size=sizes, weight=config['weight'], iter=config['iterations'], rep=config['replicates'])

rule get_dists:
	input:
		'{{input_dir}}/{dataset}.fasta'.format(dataset=dataset)
	output:
		'{input_dir}/{dataset}.dist'
	shell:
		'mothur "#set.dir(input={input_dir}, output={input_dir}); '
		'dist.seqs(fasta={dataset}.fasta, cutoff=0.03);"'

rule split_weighted_subsample:
	input:
		count="{input_dir}/{dataset}.count_table",
		dist="{input_dir}/{dataset}.dist"
	params:
		size="{size}",
		weight="{weight}",
		iter="{iter}",
		rep="{rep}"
	output:
		"{input_dir}/subsamples/size={size}_weight={weight}_i={iter}_r={rep}/sample.accnos"
	script:
		"code/analysis/weighted_subsample.R"

"""
rule optifit_test:
	shell:
		'mothur "#set.seed(seed=${SEED1}); '
		'set.dir(output=${OUTPUTDIR}, input=${INPUTDIR}); '
		'get.seqs(accnos=${OUTPUTDIR}sample.accnos, fasta=${PREFIX}${DATASET}.fasta); '
		'get.seqs(accnos=${OUTPUTDIR}sample.accnos, count=${PREFIX}${DATASET}.count_table); '
		'get.dists(column=${PREFIX}${DATASET}.dist, accnos=current); '
		'rename.file(fasta=current, count=current, accnos = current, column=current, prefix=${PREFIX}sample); '
		'set.seed(seed=${SEED2}); '
		'cluster(column=current, count=current); '
		'set.seed(seed=${SEED1}); '
		'remove.seqs(fasta=${PREFIX}${DATASET}.fasta, count=${PREFIX}${DATASET}.count_table, accnos=${OUTPUTDIR}${PREFIX}sample.accnos); '
		'list.seqs(fasta=current); '
		'get.dists(column=${PREFIX}${DATASET}.dist, accnos=current); '
		'rename.file(fasta=current, count=current, column=current, accnos=current, prefix=${PREFIX}reference); '
		'set.seed(seed=${SEED2}); '
		'cluster(column=current, count=current); '
		'set.dir(input=${OUTPUTDIR}); '
		'cluster.fit(reflist=${PREFIX}reference.opti_mcc.list, refcolumn=${PREFIX}reference.dist, refcount=${PREFIX}reference.count_table, reffasta=${PREFIX}reference.fasta, fasta=${PREFIX}sample.fasta, count=${PREFIX}sample.count_table, column=${PREFIX}sample.dist, printref=t); '
		'rename.file(file=${PREFIX}sample.optifit_mcc.sensspec, prefix=${PREFIX}sample.open.ref); '
		'cluster.fit(reflist=${PREFIX}reference.opti_mcc.list, refcolumn=${PREFIX}reference.dist, refcount=${PREFIX}reference.count_table, reffasta=${PREFIX}reference.fasta, fasta=${PREFIX}sample.fasta, count=${PREFIX}sample.count_table, column=${PREFIX}sample.dist, printref=t, method=closed); '
		'rename.file(file=${PREFIX}sample.optifit_mcc.sensspec, prefix=${PREFIX}sample.closed.ref); '
		'cluster.fit(reflist=${PREFIX}reference.opti_mcc.list, refcolumn=${PREFIX}reference.dist, refcount=${PREFIX}reference.count_table, reffasta=${PREFIX}reference.fasta, fasta=${PREFIX}sample.fasta, count=${PREFIX}sample.count_table, column=${PREFIX}sample.dist, printref=f); '
		'rename.file(file=${PREFIX}sample.optifit_mcc.sensspec, prefix=${PREFIX}sample.open.noref); '
		'cluster.fit(reflist=${PREFIX}reference.opti_mcc.list, refcolumn=${PREFIX}reference.dist, refcount=${PREFIX}reference.count_table, reffasta=${PREFIX}reference.fasta, fasta=${PREFIX}sample.fasta, count=${PREFIX}sample.count_table, column=${PREFIX}sample.dist, printref=f, method=closed); '
		'rename.file(file=${PREFIX}sample.optifit_mcc.sensspec, prefix=${PREFIX}sample.closed.noref);"'
"""
