#!/bin/bash

#Usage: optifit_test.sh INPUTDIR OUTPUTDIR DATASET SIZE WEIGHT SEED PREFIX

INPUTDIR=$1
OUTPUTDIR=$2
DATASET=$3 #Required: Dataset to use (human, mice, marine, soil)
SIZE=$4 #Required: Size is an int and tells the script how many sequences to cut out to fit against the rest of the sample
WEIGHT=$5 #Required: Boolean, whether or not to weight the subsampling by abundance
SEED1=$6 #Required: Random seed for mothur to use when cutting the data
SEED2=$7 #Required: Random seed for mothur to use when clustering
PREFIX=$8 #Optional: PREFIX allows you to add an optional PREFIX before DATASET. in case there are alternative files to use

mkdir -p ${OUTPUTDIR}

Rscript code/analysis/weighted_subsample.R ${INPUTDIR}${PREFIX}${DATASET}.count_table $SIZE $WEIGHT

mothur "#set.seed(seed=${SEED1});
	set.dir(output=${OUTPUTDIR}, input=${INPUTDIR});
	get.seqs(accnos=sample.accnos, fasta=${PREFIX}${DATASET}.fasta);
	get.seqs(accnos=sample.accnos, count=${PREFIX}${DATASET}.count_table);
	get.dists(column=${PREFIX}${DATASET}.dist, accnos=current);
	rename.file(fasta=current, count=current, accnos = current, column=current, prefix=${PREFIX}sample);
	set.seed(seed=${SEED2});
	cluster(column=current, count=current);
	set.seed(seed=${SEED1});
	remove.seqs(fasta=${PREFIX}${DATASET}.fasta, count=${PREFIX}${DATASET}.count_table, accnos=${PREFIX}sample.accnos);
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
