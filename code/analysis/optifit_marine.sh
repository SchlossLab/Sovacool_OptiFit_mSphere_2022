#!/bin/bash

#Usage: optifit_marine.sh INPUTDIR OUTPUTDIR SIZE PREFIX

INPUTDIR=$1
OUTPUTDIR=$2 #Directory to put output in
SIZE=$3 #Required: Size is an int and tells the script how many sequences to cut out to fit against the rest of the sample
PREFIX=$4 #Optional: PREFIX allows you to add an optional PREFIX after marine. in case there are alternative files to use

mkdir -p ${OUTPUTDIR}

mothur "#set.dir(output=${OUTPUTDIR});
	sub.sample(inputdir=${INPUTDIR}, fasta=${PREFIX}marine.fasta, size=$SIZE);
	list.seqs(fasta=current);
	get.seqs(accnos=current, count=${PREFIX}marine.count_table);
	get.dists(column=${PREFIX}marine.dist, accnos=current);
	rename.file(fasta=current, count=current, accnos = current, column=current, prefix=${PREFIX}sample);
	cluster(column=current, count=current);
	remove.seqs(fasta=${PREFIX}marine.fasta, count=${PREFIX}marine.count_table, accnos=current);
	list.seqs(fasta=current);
	get.dists(column=${PREFIX}marine.dist, accnos=current);
	rename.file(fasta=current, count=current, column=current, prefix=${PREFIX}reference);
	cluster(column=current, count=current);
	cluster.fit(reflist=${PREFIX}reference.opti_mcc.list, refcolumn=${PREFIX}reference.dist, refcount=${PREFIX}reference.count_table, reffasta=${PREFIX}reference.fasta, fasta=${PREFIX}sample.fasta, count=${PREFIX}sample.count_table, column=${PREFIX}sample.dist, printref=t)"
