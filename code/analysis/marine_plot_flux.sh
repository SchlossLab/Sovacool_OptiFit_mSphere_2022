#!/bin/bash

#Use this to plot the output of the flux jobs after they have finished running

#Usage: marine_pipeline.sh outputdir numseqs trimsize
OUTPUTDIR=$1 #Directory to put output in (must have trailing /)
NUMSEQS=$2 #number of seqs in the dataset (must equal trimsize if trimming)
TRIMSIZE=$3 #size to trim the original dataset too if not using the whole set

mkdir -p $OUTPUTDIR

FINAL=${OUTPUTDIR}marine.${PREFIX}sensspec.final

#Once all jobs are completed
for REFPI in {1..9}
do
	for I in {1..10} #10 iters for each REFP
	do
		REFP=$((REFPI*10)) #Counter increments by 1, but we want to increment by 10
		SEQNUM=$(($NUMSEQS-$REFP*$NUMSEQS/100)) #Calculate the actual number of sequences that will be subsampled
		LINE=$(head -2 ${OUTPUTDIR}${REFPI}_${I}/${PREFIX}sample.optifit_mcc.sensspec | tail -1) #Appends sensspec data onto a permanent file that accumulates data from all runs
		REFMCC=$(awk 'FNR==2{print $13}' ${OUTPUTDIR}${REFPI}_${I}/${PREFIX}reference.opti_mcc.sensspec) #from the second line (FNR==2) print data from the 13th column ({print $13})
		SAMPMCC=$(awk 'FNR==2{print $13}' ${OUTPUTDIR}${REFPI}_${I}/${PREFIX}sample.opti_mcc.sensspec)
		echo "$LINE	$REFP	$REFMCC	$SAMPMCC" >> $FINAL
	done
done

#Plot resulting data
Rscript code/analysis/plot_marine_sensspec.R ${OUTPUTDIR}marine.${PREFIX}sensspec.final

#Get all of the logfiles out of the main directory
mv *.logfile logfiles