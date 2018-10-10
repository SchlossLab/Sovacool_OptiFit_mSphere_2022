#!/bin/bash

#Usage: marine_pipeline.sh outputdir numseqs trimsize
OUTPUTDIR=$1 #Directory to put output in
NUMSEQS=$2 #number of seqs in the dataset that will go to optifit_multi
TRIMSIZE=$3 #size to trim the original dataset too if not using the whole thing

mkdir $OUTPUTDIR

if [ ! -f data/marine/marine.fasta ] #If raw data does not exist, get the raw data
then
	./code/data/marine.batch
fi

#Subset data (optional)
#Arg to marine_trim.sh is the number of sequences to include in subset
if [ ! -z "$2" ] #if second command line argument is not an empty string
then
	./code/data/marine_trim.sh $OUTPUTDIR $TRIMSIZE
	SUFFIX=$(echo $TRIMSIZE.)
else
	SUFFIX=""
fi

mothur "#set.dir(input=${OUTPUTDIR}, output=${OUTPUTDIR});
	dist.seqs(fasta=${SUFFIX}marine.fasta, cutoff=0.03);"

#Run optifit iteratively
./code/analysis/optifit_multi.sh $OUTPUTDIR $NUMSEQS $SUFFIX

#Plot resulting data
Rscript code/analysis/plot_marine_sensspec.R ${OUTPUTDIR}marine.${SUFFIX}sensspec.final
#REFSEQS=$(Rscript code/analysis/check_connections.R data/marine/marine.${SUFFIX}connections)

#Get all of the logfiles out of the main directory
mv *.logfile logfiles