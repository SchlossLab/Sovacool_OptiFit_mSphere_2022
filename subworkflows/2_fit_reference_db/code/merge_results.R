library(dplyr)
library(purrr)
library(readr)

merge_results <- function(key) {
  snakemake@input[[key]] %>%
    map(read_tsv) %>%
    bind_rows() %>%
    write_tsv(snakemake@output[[key]])
}

merge_results("sensspec")
merge_results("bench")
