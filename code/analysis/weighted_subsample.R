#Usage: Rscript weighted_subsample.R filename size weight
#
#Takes a count_table file and takes a random subsample based on abundance
#Used to separate a dataset into a reference and a sample for optifit
#We consider the defined subsample to be the reference, and the leftovers to be the sample
#
#Args:
# filename: path to a count_table file
# size: size of subsample to take
# weight: boolean, whether or not to weight the sample
require(dplyr)
require(readr)

args <- commandArgs(trailingOnly = TRUE)

filename <- args[1]
size <- args[2]
weight <- args[3]

#filename <- "data\\marine\\marine.count_table"
#size <- 1000
#weight <- FALSE

count_table <- readr::read_tsv(file = filename)

if (weight) {
  #Take a random subsampling weighted by total abundance
  reference <- dplyr::sample_n(count_table, size, weight = count_table$total)
  
} else {
  #Take an unweighted random sample
  reference <- dplyr::sample_n(count_table, size)
}

sample <- dplyr::filter(count_table, !(Representative_Sequence %in% reference$Representative_Sequence))

#Create .accnos files for both the sample, and the sample_complement for later use in mothur
#.accnos files contain a column of sequence names and nothing else
readr::write_tsv(dplyr::select(reference, "Representative_Sequence"), col_names = FALSE,
                 path = paste(dirname(filename), "reference.accnos", sep = "/"))
readr::write_tsv(dplyr::select(sample, "Representative_Sequence"), col_names = FALSE,
                 path = paste(dirname(filename), "sample.accnos", sep = "/"))