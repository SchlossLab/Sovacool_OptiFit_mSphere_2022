# original source: https://github.com/SchlossLab/Schloss_Cluster_PeerJ_2015/blob/master/code/uc_to_list.R
# Modified by KLS
library(tidyverse)
library(glue)
uc_to_list <- function(sorted_file_name,
                       clustered_file_name,
                       list_file_name,
                       label = 0.03) {
  
	uniqued <- read.table(file=sorted_file_name, stringsAsFactors=FALSE)

	names_first_column <- uniqued[uniqued$V1=="S", "V9"]
	names_second_column <- names_first_column

	hits <- uniqued[uniqued$V1=="H", ]

	for(i in 0:(length(names_first_column)-1)){
		dups <- paste(hits[hits$V2==i, "V9"], collapse=",")
		names_second_column[i+1] <- paste(names_second_column[i+1], dups, sep=",")
	}
	names_second_column <- gsub(",$", "", names_second_column)


	clustered <- read.table(file=clustered_file_name, stringsAsFactors=FALSE)
	clustered$sequence <- 1:nrow(clustered)

	otus <- names_second_column[clustered[clustered$V1=="S", "sequence"]]
	hits <- clustered[clustered$V1=="H", ]

	for(i in 1:nrow(hits)){
		otus[hits[i,"V2"]+1] <- paste(otus[hits[i,"V2"]+1], names_second_column[hits[i,"sequence"]], sep=",") 
	}

	#list_data <- paste(otus, collapse="\t")
	#list_data <- paste("userLabel", length(otus), list_data, sep="\t")
	num_otus <- length(otus)
	otu_ids <- sapply(seq.int(num_otus), function(x) { glue("OTU_{x}") })
	list_data <- paste(paste(c("label", "numOTUs", otu_ids),
	                             collapse = "\t"),
	                       paste(c(label, num_otus, otus),
	                             collapse = "\t"),
	                       sep = "\n")
	write(list_data, list_file_name)

}

if (exists('snakemake')) {
    sorted_filename <- snakemake@input[['sorted']]
    clustered_filename <- snakemake@input[["clustered"]]
    list_filename <- snakemake@output[["list"]]
} else {
    args <- commandArgs(trailingOnly=TRUE)
    sorted_filename <- args[1]
    clustered_filename <- args[2]
    list_filename <- args[3]
}
uc_to_list(sorted_filename, clustered_filename, list_filename)
