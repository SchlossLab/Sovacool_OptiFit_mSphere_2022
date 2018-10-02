#!/bin/bash

#Usage: marine_pipeline.sh numseqs trimsize
NUMSEQS=$1
TRIMSIZE=$2

#Get data
./code/data/marine.batch

#Subset data (optional)
#Arg to marine_trim.sh is the number of sequences to include in subet
if [ ! -z "$2" ]
then
	./code/data/marine_trim.sh $TRIMSIZE
fi


#Run optifit iteratively
./code/analysis/optifit_multi.sh $NUMSEQS

#Plot resulting data
R -e "source('code/analysis/plot_marine_sensspec.R'); make_plot()"