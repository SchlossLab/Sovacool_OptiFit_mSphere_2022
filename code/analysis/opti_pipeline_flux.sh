#!/bin/bash

#Use this to run the pipeline with each individual iteration as its own job on flux, instead of
#single back to back interations on one machine

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
if [ ! -z "${TRIMSIZE}" ] #if TRIMSIZE is not an empty string
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

if [ ! -f data/${DATASET}/${DATASET}.dist ] #If dist data doesn't exist, get dists
then
	mothur "#set.dir(input=${OUTPUTDIR}, output=${OUTPUTDIR});
dist.seqs(fasta=${PREFIX}${DATASET}.fasta, cutoff=0.03);"
fi

#Schedule one job per iteration
#Do optifit using various % of the original data as the reference
#REFP iterates from 1..20, times 5 gives us 5..100 in increments of 5
#0% is skipped because that is equivalent to just running 100% (both of which are equivalent to just running opticlust)
for REFPI in {1..20}
do
	for I in {1..10} #20 iters for each REFP
	do
		for J in {1..10}
		do
			REFP=$((REFPI*5)) #Counter increments by 1, but we want to increment by 10
			SEQNUM=$(($NUMSEQS-$REFP*$NUMSEQS/100)) #Calculate the actual number of sequences that will be subsampled
			
			#Create and fire off flux jobs
			cat optihead.pbs >> job.pbs
			echo "./code/analysis/optifit_test.sh ${OUTPUTDIR} ${OUTPUTDIR}${REFPI}_${I}_${J}/ $DATASET $SEQNUM $I $J $PREFIX" >> job.pbs #Create different output subdirectories so multiple flux jobs don't interfere with each other
			qsub -N opti_${REFPI}_${I}_${J} 	job.pbs
			rm job.pbs
		done
	done
done
