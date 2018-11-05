#!/bin/bash

#Usage: opti_pipeline.sh outputdir dataset numseqs trimsize
OUTPUTDIR=$1 #Directory to put output in (must have trailing /)
DATASET=$2 #Dataset to use (human, mice, marine, soil)
WEIGHT=$3 #Required: How to weight the subsample (none, sample_abundance, sample_dists, ref_abundance, ref_dists)
NUMSEQS=$4 #number of seqs in the dataset (must equal trimsize if trimming)
TRIMSIZE=$5 #size to trim the original dataset to if not using the whole set

mkdir -p $OUTPUTDIR

if [ ! -f data/${DATASET}/${DATASET}.fasta ] #If raw data does not exist, get the raw data
then
	./code/data/${DATASET}.batch
fi

#Subset data (optional)
#Takes the full dataset and randomly selects a subset to use in script testing
#Used to test the pipeline when doing something that takes too long on the full dataset
if [ ! -z "${TRIMSIZE}" ] #if trimsize command line argument is not an empty string
then
	mothur "#set.dir(output=${OUTPUTDIR});
	sub.sample(inputdir=data/${DATASET}/, fasta=${DATASET}.fasta, size=$NUMSEQS);
	list.seqs(fasta=current);
	get.seqs(accnos=current, count=${DATASET}.count_table);
	rename.file(accnos = current, fasta=current, count=current, prefix=$NUMSEQS.${DATASET});"
	PREFIX=$(echo $TRIMSIZE.)
else
	PREFIX=""
fi

mothur "#set.dir(input=${OUTPUTDIR}, output=${OUTPUTDIR});
	dist.seqs(fasta=${PREFIX}${DATASET}.fasta, cutoff=0.03);"

