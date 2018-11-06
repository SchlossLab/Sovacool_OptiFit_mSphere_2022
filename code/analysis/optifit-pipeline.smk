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
weight = config['weight']
sizes = [math.floor(num_seqs * i/100) for i in range(5,10,5)]
iter = range(config['iterations'])
rep = range(config['replicates'])

rule all:
	input:
		expand('{input_dir}/{dataset}.{ext}', dataset=dataset, input_dir=input_dir, ext={'fasta', 'dist', 'count_table'}),
		expand("{output_dir}/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.accnos", output_dir=output_dir, dataset=dataset, size=sizes, weight=weight, iter=iter)

rule get_dists:
	input:
		'{input_dir}/{dataset}.fasta'
	output:
		'{input_dir}/{dataset}.dist'
	shell:
		'mothur "#set.dir(input={input_dir}, output={input_dir}); '
		'dist.seqs(fasta={dataset}.fasta, cutoff=0.03);"'

rule split_weighted_subsample:
	input:
		count="{input_dir}/{{dataset}}.count_table".format(input_dir=input_dir),
		dist="{input_dir}/{{dataset}}.dist".format(input_dir=input_dir)
	params:
		size="{size}",
		weight="{weight}"
	output:
		"{output_dir}/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.accnos"
	script:
		"weighted_subsample.R"

rule prep_subsample:
	input:
		fasta="{input_dir}/{dataset}.fasta",
		accnos="{output_dir}/{dataset}_weight-{weight}_size-{size}_i-{iter}/sample.accnos"
	output:
		""
	shell:
		'mothur "#set.seed(seed={iter}}); '
		'set.dir(output={output_dir}, input={input_dir}); '
		'get.seqs(accnos={input.accnos}, fasta={input.fasta}); '
		'get.seqs(accnos={input.accnos}, count={dataset}.count_table); '
		'get.dists(column={dataset}.dist, accnos=current); '
		'rename.file(fasta=current, count=current, accnos = current, column=current, prefix=sample); '
		'remove.seqs(fasta={dataset}.fasta, count={dataset}.count_table, accnos={output_dir}sample.accnos); '
		'list.seqs(fasta=current); '
		'get.dists(column={dataset}.dist, accnos=current); '
		'rename.file(fasta=current, count=current, column=current, accnos=current, prefix=reference)"'

rule run_optifit:
	shell:
		'mothur "set.seed(seed=${rep}); '
		'cluster(column=current, count=current); '
		'set.dir(input={output_dir}); '
		'cluster.fit(reflist=reference.opti_mcc.list, refcolumn=reference.dist, refcount=reference.count_table, reffasta=reference.fasta, fasta=sample.fasta, count=sample.count_table, column=sample.dist, printref=t); '
		'rename.file(file=sample.optifit_mcc.sensspec, prefix=sample.open.ref); '
		'cluster.fit(reflist=reference.opti_mcc.list, refcolumn=reference.dist, refcount=reference.count_table, reffasta=reference.fasta, fasta=sample.fasta, count=sample.count_table, column=sample.dist, printref=t, method=closed); '
		'rename.file(file=sample.optifit_mcc.sensspec, prefix=sample.closed.ref); '
		'cluster.fit(reflist=reference.opti_mcc.list, refcolumn=reference.dist, refcount=reference.count_table, reffasta=reference.fasta, fasta=sample.fasta, count=sample.count_table, column=sample.dist, printref=f); '
		'rename.file(file=sample.optifit_mcc.sensspec, prefix=sample.open.noref); '
		'cluster.fit(reflist=reference.opti_mcc.list, refcolumn=reference.dist, refcount=reference.count_table, reffasta=reference.fasta, fasta=sample.fasta, count=sample.count_table, column=sample.dist, printref=f, method=closed); '
		'rename.file(file=sample.optifit_mcc.sensspec, prefix=sample.closed.noref)"'
