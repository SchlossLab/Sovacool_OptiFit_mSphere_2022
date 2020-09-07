library(dplyr)
library(purrr)
library(readr)

source(snakemake@input[["fcns"]])

merge_results <- function(key) {
  snakemake@input[[key]] %>%
    map(read_tsv) %>%
    bind_rows() %>%
    write_tsv(snakemake@output[[key]])
}
merge_results("tsv")
