#!/bin/bash

#Usage: marine_pipeline.sh numseqs trimsize
NUMSEQS=$1
TRIMSIZE=$2

#Get data
./code/data/marine.batch

#Subset data (optional)
#Arg to marine_trim.sh is the number of sequences to include in subet
if [ ! -z "$2" ] #if second command line argument is not an empty string
then
	./code/data/marine_trim.sh $TRIMSIZE
	SUFFIX=$(echo $TRIMSIZE.)
else
	SUFFIX=""
fi


#Run optifit iteratively
./code/analysis/optifit_multi.sh $NUMSEQS $SUFFIX

#Plot resulting data
R -e "source('code/analysis/plot_marine_sensspec.R'); make_plot()"