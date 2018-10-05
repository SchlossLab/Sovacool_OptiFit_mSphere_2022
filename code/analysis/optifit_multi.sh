#!/bin/bash

#Usage: optifit_multi.sh numseqs suffix
#numseqs is an integer argument telling the script how many total sequences are in your dataset
#suffix is an optional argument in case you are using a subset of the original data
#Takes a dataset and runs optifit on it using incremental amounts of the original dataset as a reference
#Creates a table with the sensspec data for each runs

MARINE=data/marine
FINAL=data/marine/marine.sensspec.final
NUMSEQS=$1
SUFFIX=$2

rm $FINAL
touch $FINAL
echo "iter	label	cutoff	numotus	tp	tn	fp	fn	sensitivity	specificity	ppv	npv	fdr	accuracy	mcc	f1score	refp	refmcc	sampmcc" >> $FINAL

#Do optifit using various % of the original data as the reference
#REFP iterates from 1..9, times 10 gives us 10..90 in increments of 10
#0% and 100% are skipped because that is equivalent to just running opticlust
for REFPI in {1..9}
do
	for I in {1..10} #10 iters for each REFP
	do
		REFP=$((REFPI*10)) #Counter increments by 1, but we want to increment by 10
		SEQNUM=$(($NUMSEQS-$REFP*$NUMSEQS/100)) #Calculate the actual number of sequences that will be subsampled
		./code/analysis/optifit_marine.sh $SEQNUM $SUFFIX
		LINE=$(head -2 $MARINE/sample.optifit_mcc.sensspec | tail -1) #Appends sensspec data onto a permanent file that accumulates data from all runs
		#REFSEQS=$(Rscript code/analysis/check_connections.R data/marine/marine.${SUFFIX}connections)
		REFMCC=$(awk 'FNR==2{print $13}' data/marine/reference.opti_mcc.sensspec) #from the second line (FNR==2) print data from the 13th column ({print $13})
		SAMPMCC=$(awk 'FNR==2{print $13}' data/marine/sample.opti_mcc.sensspec)
		echo "$LINE	$REFP	$REFMCC	$SAMPMCC" >> $FINAL
	done
done
