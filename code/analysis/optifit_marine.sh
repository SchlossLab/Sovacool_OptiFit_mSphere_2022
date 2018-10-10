#!/bin/bash

#Usage: optifit_marine.sh size suffix

OUTPUTDIR=$1 #Directory to put output in
SIZE=$2 #Required: Size is an int and tells the script how many sequences to cut out to fit against the rest of the sample
SUFFIX=$3 #Optional: Suffix allows you to add an optional suffix after marine. in case there are alternative files to use

mothur "#set.dir(output=${OUTPUTDIR});
	sub.sample(inputdir=${OUTPUTDIR}, fasta=${SUFFIX}marine.fasta, size=$SIZE);
	list.seqs(fasta=current);
	get.seqs(accnos=current, count=${SUFFIX}marine.count_table);
	get.dists(column=${SUFFIX}marine.dist, accnos=current);
	rename.file(fasta=current, count=current, accnos = current, column=current, prefix=${SUFFIX}sample);
	cluster(column=current, count=current);
	remove.seqs(fasta=${SUFFIX}marine.fasta, count=${SUFFIX}marine.count_table, accnos=current);
	list.seqs(fasta=current);
	get.dists(column=${SUFFIX}marine.dist, accnos=current);
	rename.file(fasta=current, count=current, column=current, prefix=${SUFFIX}reference);
	cluster(column=current, count=current);
	cluster.fit(reflist=${SUFFIX}reference.opti_mcc.list, refcolumn=${SUFFIX}reference.dist, refcount=${SUFFIX}reference.count_table, reffasta=${SUFFIX}reference.fasta, fasta=${SUFFIX}sample.fasta, count=${SUFFIX}sample.count_table, column=${SUFFIX}sample.dist, printref=t)"
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	