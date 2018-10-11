#!/bin/bash

#Use this to plot the output of the flux jobs after they have finished running

#Usage: marine_pipeline.sh OUTPUTDIR PREFIX
OUTPUTDIR=$1 #Directory to put output in (must have trailing /)
PREFIX=$2 #Prefix used to create sensspec files

mkdir -p $OUTPUTDIR
FINAL=${OUTPUTDIR}marine.${PREFIX}sensspec.final

touch $FINAL
echo "iter	label	cutoff	numotus	tp	tn	fp	fn	sensitivity	specificity	ppv	npv	fdr	accuracy	mcc	f1score	refp	refmcc	sampmcc" >> $FINAL

#Once all jobs are completed
for REFPI in {1..9}
do
	for I in {1..10} #10 iters for each REFP
	do
		REFP=$((REFPI*10)) #Counter increments by 1, but we want to increment by 10
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