#!/bin/bash

#Usage: marine_trim.sh numseqs
#numseqs is a required integer argument that tells the script how many sequences to take in the subsample
#Takes the full marine dataset and randomly selects a subset to use in script testing
#Original dataset is too big so stuff takes forever to run on it
#Only using the subset if data to test pipeline, use full dataset for actual analysis

#Requires Mothur installed in path

MARINE=data/marine
NUMSEQS=$1

mothur "#set.dir(output=$MARINE);
	sub.sample(inputdir=$MARINE, fasta=marine.fasta, size=$NUMSEQS);
	list.seqs(fasta=current);
	get.seqs(accnos=current, count=marine.count_table);
	rename.file(fasta=current, count=current, prefix=marine.$NUMSEQS);"
