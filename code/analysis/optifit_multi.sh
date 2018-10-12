#!/bin/bash

#Usage: optifit_multi.sh outputdir numseqs PREFIX
OUTPUTDIR=$1 #Directory to put output in
NUMSEQS=$2 #numseqs is an integer argument telling the script how many total sequences are in your dataset
PREFIX=$3 #PREFIX is an optional argument in case you are using a subset of the original data

#Takes a dataset and runs optifit on it using incremental amounts of the original dataset as a reference
#Creates a table with the sensspec data for each runs

FINAL=${OUTPUTDIR}marine.${PREFIX}sensspec.final


rm $FINAL
touch $FINAL
echo "iter	label	cutoff	numotus	tp	tn	fp	fn	sensitivity	specificity	ppv	npv	fdr	accuracy	mcc	f1score	refp	refpi	type" >> $FINAL

#Do optifit using various % of the original data as the reference
#REFP iterates from 1..9, times 10 gives us 10..90 in increments of 10
#0% and 100% are skipped because that is equivalent to just running opticlust
for REFPI in {1..2}
do
	for I in {1..2} #10 iters for each REFP
	do
		REFP=$((REFPI*10)) #Counter increments by 1, but we want to increment by 10
		SEQNUM=$(($NUMSEQS-$REFP*$NUMSEQS/100)) #Calculate the actual number of sequences that will be subsampled
		./code/analysis/optifit_marine.sh $OUTPUTDIR $OUTPUTDIR $SEQNUM $PREFIX
		
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
