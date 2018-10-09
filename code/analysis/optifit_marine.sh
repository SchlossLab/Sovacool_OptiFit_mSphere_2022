#!/bin/bash

#Usage: optifit_marine.sh size suffix

SIZE=$1 #Required: Size is an int and tells the script how many sequences to cut out to fit against the rest of the sample
SUFFIX=$2 #Optional: Suffix allows you to add an optional suffix after marine. in case there are alternative files to use

MARINE=data/marine

mothur "#set.dir(output=$MARINE);
	sub.sample(inputdir=$MARINE, fasta=marine.${SUFFIX}fasta, size=$SIZE);
	list.seqs(fasta=current);
	get.seqs(accnos=current, count=marine.${SUFFIX}count_table);
	get.dists(column=marine.${SUFFIX}dist, accnos=current);
	rename.file(fasta=current, count=current, accnos = current, column=current, prefix=${SUFFIX}sample);
	cluster(column=current, count=current);
	remove.seqs(fasta=marine.${SUFFIX}fasta, count=marine.${SUFFIX}count_table, accnos=current);
	list.seqs(fasta=current);
	get.dists(column=marine.${SUFFIX}dist, accnos=current);
	rename.file(fasta=current, count=current, column=current, prefix=${SUFFIX}reference);
	cluster(column=current, count=current);
	cluster.fit(reflist=${SUFFIX}reference.opti_mcc.list, refcolumn=${SUFFIX}reference.dist, refcount=${SUFFIX}reference.count_table, reffasta=${SUFFIX}reference.fasta, fasta=${SUFFIX}sample.fasta, count=${SUFFIX}sample.count_table, column=${SUFFIX}sample.dist, inputdir=$MARINE, printref=t)"