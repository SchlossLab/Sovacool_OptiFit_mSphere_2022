#!/bin/bash

#Use this to plot the output of the flux jobs after they have finished running

#Usage: marine_pipeline.sh OUTPUTDIR PREFIX
OUTPUTDIR=$1 #Directory to put output in (must have trailing /)
PREFIX=$2 #Prefix used to create sensspec files

mkdir -p $OUTPUTDIR
FINAL=${OUTPUTDIR}${PREFIX}marine.sensspec.final

touch $FINAL
echo "iter	label	cutoff	numotus	tp	tn	fp	fn	sensitivity	specificity	ppv	npv	fdr	accuracy	mcc	f1score	refp	refpi	type" >> $FINAL

#Once all jobs are completed
for REFPI in {1..9}
do
	for I in {1..10} #10 iters for each REFP
	do
		REFP=$((REFPI*10)) #Counter increments by 1, but we want to increment by 10
		#opticlust_marine.sh will run opticlust on the reference alone and the sample alone, and then
		#optifit with fitting the sample to the reference with all pairwise possibilities of
		#method = (open, closed) and printref = (T, F)
		
		#REF and SAMP are run with opticlust, which has a output with a different format than optifit output
		#Optifit output is missing the iter and numOTUs columns, so we use sed {'s_\(\S*\t\S*\t\)\(.*\)_\1\t\2_'} to insert a tab in the beginning of the line
		#and then find the second tab and turn it into two tabs to give us the correct number of columns
		#head -2: return first 2 lines, tail -1 from those 2 lines take only the bottom one
		REF=$(head -2 ${OUTPUTDIR}${PREFIX}reference.opti_mcc.sensspec | tail -1 | sed 's_\(\S*\t\S*\t\)\(.*\)_\t\1\t\2_')
		SAMP=$(head -2 ${OUTPUTDIR}${PREFIX}sample.opti_mcc.sensspec | tail -1 | sed 's_\(\S*\t\S*\t\)\(.*\)_\t\1\t\2_')
		SAMP_O_REF=$(head -2 ${OUTPUTDIR}${PREFIX}sample.open.ref.sensspec | tail -1)
		SAMP_C_REF=$(head -2 ${OUTPUTDIR}${PREFIX}sample.closed.ref.sensspec | tail -1)
		
		#REF and SAMP were run with opticlust, which produces sensspec files with 2 less columns than optifit
		#Add two extra tabs at the beginning of their lines so that confusion matrix values line up
		#REF and SAMP also have records that end in a tab, so one less tab at the end
		echo "${REF}${REFP}	$I	REF" >> $FINAL
		echo "${SAMP}${REFP}	$I	SAMP" >> $FINAL
		echo "$SAMP_O_REF	$REFP	$I	SAMP_O_REF" >> $FINAL
		echo "$SAMP_C_REF	$REFP	$I	SAMP_C_REF" >> $FINAL
	done
done

#Plot resulting data
Rscript code/analysis/plot_marine_sensspec.R ${OUTPUTDIR}marine.${PREFIX}sensspec.final

#Get all of the logfiles out of the main directory
mv *.logfile logfiles