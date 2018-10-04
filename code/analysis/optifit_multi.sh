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
echo "iter	label	cutoff	numotus	tp	tn	fp	fn	sensitivity	specificity	ppv	npv	fdr	accuracy	mcc	f1score	refp" >> $FINAL

#Do optifit using various % of the original data as the reference
#REFP iterates from 1..19, times 5 gives us 5..95 in increments of 5
#0% and 100% are skipped because that is equivalent to just running opticlust
for REFPI in {1..19}
do
	for I in {1..10} #10 iters for each REFP
	do
		REFP=$((REFPI*5))
		#Calculate the actual number of sequences that will be subsampled
		SEQNUM=$(($NUMSEQS-$REFP*$NUMSEQS/100))
		./code/analysis/optifit_marine.sh $SEQNUM $SUFFIX
		#Appends sensspec data onto a permanent file that accumulates data from all runs
		LINE=$(head -2 $MARINE/sample.optifit_mcc.sensspec | tail -1)
		echo "$LINE	$REFP" >> $FINAL
		#Adds the percent used as reference to the last column of the last line of the file
		#sed '$s/$/,' "\t$REFP" '/' $FINAL | tee $MARINE/temp.txt
		#mv $MARINE/temp.txt $FINAL
	done
done
