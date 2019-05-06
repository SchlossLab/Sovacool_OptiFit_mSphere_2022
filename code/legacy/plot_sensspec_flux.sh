#!/bin/bash

#Use this to plot the output of the flux jobs after they have finished running

#Usage: sensspec_plot_flux.sh OUTPUTDIR DATASET PREFIX
OUTPUTDIR=$1 #Directory where opti_pipeline_flux output to (must have trailing /)
DATASET=$2 #Dataset to use (human, mice, marine, soil)
PREFIX=$3 #Prefix used to create sensspec files

mkdir -p $OUTPUTDIR
FINAL=${OUTPUTDIR}${PREFIX}${DATASET}.sensspec.final

rm $FINAL
touch $FINAL
echo "iter	label	cutoff	numotus	tp	tn	fp	fn	sensitivity	specificity	ppv	npv	fdr	accuracy	mcc	f1score	refp	refpi	refpj	type" >> $FINAL

#Once all jobs are completed
for REFPI in {1..20} #Percent of data to use as reference
do
	for I in {1..10} #Number of ways to cut the data at each ref percent
	do
		for J in {1..10} #Number of times to run opticlust/optifit per cut
		do
			REFP=$((REFPI*5)) #Counter increments by 1, but we want to increment by 5
			#optifit_test.sh will run opticlust on the reference alone and the sample alone, and then
			#optifit with fitting the sample to the reference with all pairwise possibilities of
			#method = (open, closed) and printref = (T, F)
			
			#REF and SAMP are run with opticlust, which has a output with a different format than optifit output
			#Optifit output is missing the iter and numOTUs columns, so we use sed {'s_\(\S*\t\S*\t\)\(.*\)_\1\t\2_'} to insert a tab in the beginning of the line
			#and then find the second tab and turn it into two tabs to give us the correct number of columns
			#head -2: return first 2 lines, tail -1 from those 2 lines take only the bottom one
			REF=$(head -2 ${OUTPUTDIR}${REFPI}_${I}_${J}/${PREFIX}reference.opti_mcc.sensspec | tail -1 | sed 's_\(\S*\t\S*\t\)\(.*\)_\t\1\t\2_')
			SAMP=$(head -2 ${OUTPUTDIR}${REFPI}_${I}_${J}/${PREFIX}sample.opti_mcc.sensspec | tail -1 | sed 's_\(\S*\t\S*\t\)\(.*\)_\t\1\t\2_')
			SAMP_O_REF=$(head -2 ${OUTPUTDIR}${REFPI}_${I}_${J}/${PREFIX}sample.open.ref.sensspec | tail -1)
			SAMP_C_REF=$(head -2 ${OUTPUTDIR}${REFPI}_${I}_${J}/${PREFIX}sample.closed.ref.sensspec | tail -1)
			SAMP_O_NOREF=$(head -2 ${OUTPUTDIR}${REFPI}_${I}_${J}/${PREFIX}sample.open.noref.sensspec | tail -1)
			SAMP_C_NOREF=$(head -2 ${OUTPUTDIR}${REFPI}_${I}_${J}/${PREFIX}sample.closed.noref.sensspec | tail -1)
			
			#REF and SAMP were run with opticlust, which produces sensspec files with 2 less columns than optifit
			#Add two extra tabs at the beginning of their lines so that confusion matrix values line up
			#REF and SAMP also have records that end in a tab, so one less tab at the end
			echo "${REF}${REFP}	$I	$J	REF" >> $FINAL
			echo "${SAMP}${REFP}	$I	$J	SAMP" >> $FINAL
			echo "$SAMP_O_REF	$REFP	$I	$J	SAMP_O_REF" >> $FINAL
			echo "$SAMP_C_REF	$REFP	$I	$J	SAMP_C_REF" >> $FINAL
			echo "$SAMP_O_NOREF	$REFP	$I	$J	SAMP_O_NOREF" >> $FINAL
			echo "$SAMP_C_NOREF	$REFP	$I	$J	SAMP_C_NOREF" >> $FINAL
		done
	done
done

#Plot resulting data
Rscript code/analysis/plot_sensspec.R ${OUTPUTDIR}${DATASET}.${PREFIX}sensspec.final ${DATASET}

#Get all of the logfiles out of the main directory
mv *.logfile logfiles