#!/bin/bash

#Usage: marine_pipeline.sh outputdir numseqs trimsize
OUTPUTDIR=$1 #Directory to put output in (must have trailing /)
DATASET=$2 #Dataset to use (human, mice, marine, soil)
NUMSEQS=$3 #number of seqs in the dataset (must equal trimsize if trimming)
TRIMSIZE=$4 #size to trim the original dataset too if not using the whole set

mkdir -p $OUTPUTDIR

if [ ! -f data/${DATASET}/${DATASET}.fasta ] #If raw data does not exist, get the raw data
then
	./code/data/${DATASET}.batch
fi

#Subset data (optional)
#Takes the full dataset and randomly selects a subset to use in script testing
#Used to test the pipeline when doing something that takes too long on the full dataset
if [ ! -z "${TRIMSIZE}" ] #if second command line argument is not an empty string
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

#Run optifit iteratively
./code/analysis/optifit_multi.sh $OUTPUTDIR $DATASET $NUMSEQS $PREFIX

#Plot resulting data
Rscript code/analysis/plot_sensspec.R ${OUTPUTDIR}${PREFIX}${DATASET}.sensspec.final ${DATASET}

#Get all of the logfiles out of the main directory
mv *.logfile logfiles