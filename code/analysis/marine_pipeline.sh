#!/bin/bash

#Usage: marine_pipeline.sh outputdir numseqs trimsize
OUTPUTDIR=$1 #Directory to put output in (must have trailing /)
NUMSEQS=$2 #number of seqs in the dataset (must equal trimsize if trimming)
TRIMSIZE=$3 #size to trim the original dataset too if not using the whole set

mkdir -p $OUTPUTDIR

if [ ! -f data/marine/marine.fasta ] #If raw data does not exist, get the raw data
then
	./code/data/marine.batch
fi

#Subset data (optional)
#Arg to marine_trim.sh is the number of sequences to include in subset
if [ ! -z "$2" ] #if second command line argument is not an empty string
then
	./code/data/marine_trim.sh $OUTPUTDIR $TRIMSIZE
	PREFIX=$(echo $TRIMSIZE.)
else
	PREFIX=""
fi

mothur "#set.dir(input=${OUTPUTDIR}, output=${OUTPUTDIR});
	dist.seqs(fasta=${PREFIX}marine.fasta, cutoff=0.03);"

#Run optifit iteratively
./code/analysis/optifit_multi.sh $OUTPUTDIR $NUMSEQS $PREFIX

#Plot resulting data
Rscript code/analysis/plot_marine_sensspec.R ${OUTPUTDIR}marine.${PREFIX}sensspec.final

#Get all of the logfiles out of the main directory
mv *.logfile logfiles