#!/bin/bash

#Usage: optifit_marine.sh size suffix
#Required: Size is an int and tells the script how many sequences to cut out to fit against the rest of the sample
#Optional: Suffix allows you to add an optional suffix after marine. in case there are alternative files to use

MARINE=data/marine

#select a small sample to fit
#sub.sample(inputdir=$MARINE, fasta=marine.$2fasta, size=$1)
#whose in the sample
#list.seqs(fasta=current)
#select sample from count file
#get.seqs(accnos=current, count=marine.$2count_table);
#rename small sample
#rename.file(fasta=current, count=current, accnos = current, prefix=sample);
#get dists for sample
#dist.seqs(fasta=current, cutoff=0.03)
#create reference by removing sample
#remove.seqs(fasta=marine.$2fasta, count=marine.$2count_table, accnos=current);
#rename ref
#rename.file(fasta=current, count=current, prefix=reference)
#get dists for reference
#dist.seqs(fasta=current, cutoff=0.03)
#cluster reference
#cluster(column=current, count=current)
#fit sample to reference
#cluster.fit(reflist=reference.opti_mcc.list, refcolumn=reference.dist, refcount=reference.count_table, reffasta=reference.fasta, fasta=sample.fasta, count=sample.count_table, column=sample.dist, inputdir=$MARINE, printref=t)"

mothur "#set.dir(output=$MARINE);
	sub.sample(inputdir=$MARINE, fasta=marine.$2fasta, size=$1);
	list.seqs(fasta=current);
	get.seqs(accnos=current, count=marine.$2count_table);
	rename.file(fasta=current, count=current, accnos = current, prefix=sample);
	dist.seqs(fasta=current, cutoff=0.03);
	remove.seqs(fasta=marine.$2fasta, count=marine.$2count_table, accnos=current);
	rename.file(fasta=current, count=current, prefix=reference);
	dist.seqs(fasta=current, cutoff=0.03);
	cluster(column=current, count=current);
	cluster.fit(reflist=reference.opti_mcc.list, refcolumn=reference.dist, refcount=reference.count_table, reffasta=reference.fasta, fasta=sample.fasta, count=sample.count_table, column=sample.dist, inputdir=$MARINE, printref=t)"