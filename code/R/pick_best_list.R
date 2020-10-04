library(dplyr)
library(readr)
source(snakemake@input[["fcns"]])
log_smk()

sensspec <- snakemake@input[["sensspec"]]
best_iter <- sensspec %>%
  read_tsv() %>%
  dplyr::top_n(1, mcc) %>%
  pull(iter) %>%
  .[1]
best_list_file <- snakemake@input[["lists"]][[best_iter]]
message(paste("Copying best list file:", best_list_file))
file.copy(best_list_file, snakemake@output[["list"]])
