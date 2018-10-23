#!/bin/bash

#Usage: trim.sh outputdir numseqs
OUTPUTDIR=$1 #Directory to put output in
DATASET=$2 #Dataset to use (human, mice, marine, soil)
NUMSEQS=$3 #numseqs is a required integer argument that tells the script how many sequences to take in the subsample

#Takes the full dataset and randomly selects a subset to use in script testing
#Original dataset is too big so stuff takes forever to run on it
#Only using the subset if data to test pipeline, use full dataset for actual analysis

#Requires Mothur installed in path

mothur "#set.dir(output=${OUTPUTDIR});
	sub.sample(inputdir=data/${DATASET}/, fasta=${DATASET}.fasta, size=$NUMSEQS);
	list.seqs(fasta=current);
	get.seqs(accnos=current, count=${DATASET}.count_table);
	rename.file(accnos = current, fasta=current, count=current, prefix=$NUMSEQS.${DATASET});"
