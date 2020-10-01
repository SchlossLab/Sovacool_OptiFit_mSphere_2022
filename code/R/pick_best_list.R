library(dplyr)
library(readr)
source(snakemake@input[['fcns']])
smk_log()

sensspec <- snakemake@input[['sensspec']]
max_iter <- sensspec %>% 
  read_tsv() %>% 
  slice_max(order_by = mcc) %>% 
  pull(iter)
best_list_file <- snakemake@input[['lists']][[max_iter]]
message(paste('Copying best list file:', best_list_file))
file.copy(best_list_file, snakemake@output[['list']])