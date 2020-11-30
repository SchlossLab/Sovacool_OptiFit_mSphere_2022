# original source: https://github.com/SchlossLab/Schloss_Cluster_PeerJ_2015/blob/master/code/uc_to_list.R
# Modified Nov. 2020 by KLS
library(tidyverse)
uc_to_list <- function(clustered_file_name,
                       list_file_name,
                       label = 0.03) {
  clustered <- read.table(
    file = clustered_file_name,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      seq_id = str_replace(V9, "(.*);size=\\d+$", "\\1"),
      otu_id = paste0("OTU_", V2),
      record_type = V1
    ) %>%
    filter(record_type %in% c("S", "H")) %>%
    select(seq_id, otu_id)
  otu_ids <- clustered %>%
    pull(otu_id) %>%
    unique()
  num_otus <- length(otu_ids)
  otus <- sapply(otu_ids, function(x) {
    clustered %>%
      filter(otu_id == x) %>%
      pull(seq_id) %>%
      paste(collapse = ",")
  })
  # storing this in memory as a giant string is a bad idea,
  # but I only have to run it once per vsearch clustering job,
  # so I don't really care.
  list_data <- paste(paste(c("label", "numOTUs", otu_ids),
    collapse = "\t"
  ),
  paste(c(label, num_otus, otus),
    collapse = "\t"
  ),
  sep = "\n"
  )
  write(list_data, list_file_name)
}

uc_to_list(snakemake@input[["clustered"]], snakemake@output[["list"]])
