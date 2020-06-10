library(dplyr)
library(readr)

mutate_columns <- function(df) {
  df %>% mutate(
    dataset = snakemake@params[["dataset"]],
    ref = NA,
    region = NA,
    seed = snakemake@params[["seed"]],
    method = "de_novo",
    printref = NA,
    iter = NA,
    numotus = NA
  )
}

reformat <- function(key) {
  read_tsv(snakemake@input[[key]]) %>%
    mutate_columns() %>%
    write_tsv(snakemake@output[[key]])
}

reformat("sensspec")
reformat("bench")
