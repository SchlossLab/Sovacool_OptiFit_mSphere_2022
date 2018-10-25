#!/bin/bash

#Usage: optifit_test.sh INPUTDIR OUTPUTDIR DATASET SIZE SEED PREFIX

INPUTDIR=$1
OUTPUTDIR=$2 #Directory to put output in
DATASET=$3
SIZE=$4 #Required: Size is an int and tells the script how many sequences to cut out to fit against the rest of the sample
SEED1=$5 #Required: Random seed for mothur to use when cutting the data
SEED2=$6 #Required: Random seed for mothur to use when clustering
PREFIX=$7 #Optional: PREFIX allows you to add an optional PREFIX before DATASET. in case there are alternative files to use

mkdir -p ${OUTPUTDIR}

mothur "#set.seed(seed=${SEED1});
	set.dir(output=${OUTPUTDIR});
	sub.sample(inputdir=${INPUTDIR}, fasta=${PREFIX}${DATASET}.fasta, size=$SIZE);
	list.seqs(fasta=current);
	get.seqs(accnos=current, count=${PREFIX}${DATASET}.count_table);
	get.dists(column=${PREFIX}${DATASET}.dist, accnos=current);
	rename.file(fasta=current, count=current, accnos = current, column=current, prefix=${PREFIX}sample);
	set.seed(seed=${SEED2});
	cluster(column=current, count=current);
	set.seed(seed=${SEED1});
	remove.seqs(fasta=${PREFIX}${DATASET}.fasta, count=${PREFIX}${DATASET}.count_table, accnos=current);
	list.seqs(fasta=current);
	get.dists(column=${PREFIX}${DATASET}.dist, accnos=current);
	rename.file(fasta=current, count=current, column=current, accnos=current, prefix=${PREFIX}reference);
	set.seed(seed=${SEED2});
	cluster(column=current, count=current);
	cluster.fit(reflist=${PREFIX}reference.opti_mcc.list, refcolumn=${PREFIX}reference.dist, refcount=${PREFIX}reference.count_table, reffasta=${PREFIX}reference.fasta, fasta=${PREFIX}sample.fasta, count=${PREFIX}sample.count_table, column=${PREFIX}sample.dist, printref=t);
	rename.file(file=${PREFIX}sample.optifit_mcc.sensspec, prefix=${PREFIX}sample.open.ref);
	cluster.fit(reflist=${PREFIX}reference.opti_mcc.list, refcolumn=${PREFIX}reference.dist, refcount=${PREFIX}reference.count_table, reffasta=${PREFIX}reference.fasta, fasta=${PREFIX}sample.fasta, count=${PREFIX}sample.count_table, column=${PREFIX}sample.dist, printref=t, method=closed);
	rename.file(file=${PREFIX}sample.optifit_mcc.sensspec, prefix=${PREFIX}sample.closed.ref);
	cluster.fit(reflist=${PREFIX}reference.opti_mcc.list, refcolumn=${PREFIX}reference.dist, refcount=${PREFIX}reference.count_table, reffasta=${PREFIX}reference.fasta, fasta=${PREFIX}sample.fasta, count=${PREFIX}sample.count_table, column=${PREFIX}sample.dist, printref=f);
	rename.file(file=${PREFIX}sample.optifit_mcc.sensspec, prefix=${PREFIX}sample.open.noref);
	cluster.fit(reflist=${PREFIX}reference.opti_mcc.list, refcolumn=${PREFIX}reference.dist, refcount=${PREFIX}reference.count_table, reffasta=${PREFIX}reference.fasta, fasta=${PREFIX}sample.fasta, count=${PREFIX}sample.count_table, column=${PREFIX}sample.dist, printref=f, method=closed);
	rename.file(file=${PREFIX}sample.optifit_mcc.sensspec, prefix=${PREFIX}sample.closed.noref);"
