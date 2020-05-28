library(dplyr)
library(purrr)
library(readr)

merge_results <- function(key) {
  infilename <- snakemake@input[[key]]
  dfs <- map(infilenames, read_tsv)
  bind_rows(dfs) %>% write_tsv(snakemake@output[[key]])
}
merge_results('sensspec')
merge_results('bench')