#Usage: call from snakemake script directive
#Takes a count_table file and takes a random subsample based on abundance
#Used to separate a dataset into a reference and a sample for optifit
#
require(dplyr)
require(readr)

filename_counts <- snakemake@input[['count']]
filename_dists <- snakemake@input[['dist']]
filename_output <- snakemake@output[[1]]
size_fraction <- as.numeric(snakemake@params[['size']])
weight <- snakemake@params[['weight']]

#filename <- "data\\mice\\mice.count_table"
#size <- 10000
#weight <- FALSE
#dist_file <- "data\\mice\\mice.dist"

count_table <- readr::read_tsv(file = filename_counts)

if (weight == "none") {
  #Take an unweighted random sample
  sample <- dplyr::sample_frac(count_table, size_fraction)
} else if (weight == "sample-abundance") {
  #Take a random subsampling weighted by total abundance
  sample <- dplyr::sample_frac(count_table, size_fraction, weight = count_table$total)
} else if (weight == "ref-abundance") {
  #Take a random subsampling of the complement of the sample by total abundance, and use the leftovers as sample
  sample_complement <- dplyr::sample_frac(count_table, size = 1-size_fraction, weight = count_table$total)
  sample <- dplyr::filter(count_table, !(Representative_Sequence %in% sample_complement$Representative_Sequence))
} else if (weight == "sample-dists") {
  #Take a random subsample weighted by number of pairwise connections to other seqs
  dists <- readr::read_delim(file = filename_dists, delim = " ", col_names = FALSE)

  #Each row of a dist file is a pair of sequences that are pairwise close under our threshold
  #Since we want the total number of distances for each sequence, we need to check both columns
  #for a given sequence. Do this by simply stacking both columns and counting the number of occurences
  #of each sequence. We also stack in count_table$Representative_Sequence so that every sequence appears
  #at least once, otherwise many would have a weight of zero. This means that if every sequence has X number
  #of connections, we will assign it a weight of X+1 when subsampling.
  dist_seqs <- tibble(Representative_Sequence = c(dists$X1, dists$X2, count_table$Representative_Sequence)) %>%
    group_by(Representative_Sequence) %>% #Group so that mutate knows what to count
    mutate(count = n()) %>%
    unique() %>% #Don't need duplicate rows
    ungroup() #Have to remove groups to do sample_n

  sample <- dplyr::sample_frac(dist_seqs, size_fraction, weight = dist_seqs$count)
} else if(weight == "ref-dists") {
  #Take a random subsample of the complement of the sample weighted by number of pairwise connections to other seqs,
  #and use the leftovers as the sample
  dists <- readr::read_delim(file = filename_dists, delim = " ", col_names = FALSE)

  dist_seqs <- tibble(Representative_Sequence = c(dists$X1, dists$X2, count_table$Representative_Sequence)) %>%
    group_by(Representative_Sequence) %>% #Group so that mutate knows what to count
    mutate(count = n()) %>%
    unique() %>% #Don't need duplicate rows
    ungroup() #Have to remove groups to do sample_n

  sample_complement <- dplyr::sample_frac(dist_seqs, size = 1-size_fraction, weight = dist_seqs$count)
  sample <- dplyr::filter(dist_seqs, !(Representative_Sequence %in% sample_complement$Representative_Sequence))
}

reference <- dplyr::filter(count_table, !(Representative_Sequence %in% sample$Representative_Sequence))

#Create .accnos files for both the sample, and the sample_complement for later use in mothur
#.accnos files contain a column of sequence names and nothing else
readr::write_tsv(dplyr::select(sample, "Representative_Sequence"), col_names = FALSE,
                 path = paste(filename_output, sep = "/"))
