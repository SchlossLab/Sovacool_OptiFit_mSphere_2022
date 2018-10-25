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
#Arg to marine_trim.sh is the number of sequences to include in subset
if [ ! -z "${TRIMSIZE}" ] #if second command line argument is not an empty string
then
	./code/data/trim.sh $OUTPUTDIR $DATASET $TRIMSIZE
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