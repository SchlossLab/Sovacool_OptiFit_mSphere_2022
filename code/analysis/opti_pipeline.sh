#!/bin/bash

#Usage: opti_pipeline.sh outputdir dataset numseqs trimsize
OUTPUTDIR=$1 #Directory to put output in (must have trailing /)
DATASET=$2 #Dataset to use (human, mice, marine, soil)
WEIGHT=$3 #Boolean, whether or not to weight the reference subsample
NUMSEQS=$4 #number of seqs in the dataset (must equal trimsize if trimming)
TRIMSIZE=$5 #size to trim the original dataset to if not using the whole set

mkdir -p $OUTPUTDIR

if [ ! -f data/${DATASET}/${DATASET}.fasta ] #If raw data does not exist, get the raw data
then
	./code/data/${DATASET}.batch
fi

#Subset data (optional)
#Takes the full dataset and randomly selects a subset to use in script testing
#Used to test the pipeline when doing something that takes too long on the full dataset
if [ ! -z "${TRIMSIZE}" ] #if trimsize command line argument is not an empty string
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

mothur "#set.dir(input=${OUTPUTDIR}, output=${OUTPUTDIR});
	dist.seqs(fasta=${PREFIX}${DATASET}.fasta, cutoff=0.03);"

#Run optifit iteratively
#Takes a dataset and runs optifit on it using incremental amounts of the original dataset as a reference
#Creates a table with the sensspec data for each runs

FINAL=${OUTPUTDIR}${PREFIX}${DATASET}.sensspec.final

rm $FINAL
touch $FINAL
echo "iter	label	cutoff	numotus	tp	tn	fp	fn	sensitivity	specificity	ppv	npv	fdr	accuracy	mcc	f1score	refp	refpi	iter	type" >> $FINAL

#Do optifit using various % of the original data as the reference
#REFP iterates from 1..9, times 10 gives us 10..90 in increments of 10
#0% and 100% are skipped because that is equivalent to just running opticlust
for REFPI in {5..5}
do
	for I in {1..3} #10 iters for each REFP
	do
		for J in {1..5}
		do
			REFP=$((REFPI*10)) #Counter increments by 1, but we want to increment by 10
			SEQNUM=$(($NUMSEQS-$REFP*$NUMSEQS/100)) #Calculate the actual number of sequences that will be subsampled
			./code/analysis/optifit_test.sh $OUTPUTDIR $OUTPUTDIR $DATASET $SEQNUM $WEIGHT $I $J $PREFIX
			
			#optifit_test.sh will run opticlust on the reference alone and the sample alone, and then
			#optifit with fitting the sample to the reference with all pairwise possibilities of
			#method = (open, closed) and printref = (T, F)
			
			#REF and SAMP are run with opticlust, which has an output with a different format than optifit output
			#Opticlust output is missing the iter and numOTUs columns, so we use sed {'s_\(\S*\t\S*\t\)\(.*\)_\1\t\2_'} to insert a tab in the beginning of the line
			#and then find the second tab and turn it into two tabs to give us the correct number of columns
			#head -2: return first 2 lines, tail -1 from those 2 lines take only the bottom one
			REF=$(head -2 ${OUTPUTDIR}${PREFIX}reference.opti_mcc.sensspec | tail -1 | sed 's_\(\S*\t\S*\t\)\(.*\)_\t\1\t\2_')
			SAMP=$(head -2 ${OUTPUTDIR}${PREFIX}sample.opti_mcc.sensspec | tail -1 | sed 's_\(\S*\t\S*\t\)\(.*\)_\t\1\t\2_')
			SAMP_O_REF=$(head -2 ${OUTPUTDIR}${PREFIX}sample.open.ref.sensspec | tail -1)
			SAMP_C_REF=$(head -2 ${OUTPUTDIR}${PREFIX}sample.closed.ref.sensspec | tail -1)
			SAMP_O_NOREF=$(head -2 ${OUTPUTDIR}${PREFIX}sample.open.noref.sensspec | tail -1)
			SAMP_C_NOREF=$(head -2 ${OUTPUTDIR}${PREFIX}sample.closed.noref.sensspec | tail -1)
			
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
Rscript code/analysis/plot_sensspec.R ${OUTPUTDIR}${PREFIX}${DATASET}.sensspec.final ${DATASET}

#Get all of the logfiles out of the main directory
mv *.logfile logfiles