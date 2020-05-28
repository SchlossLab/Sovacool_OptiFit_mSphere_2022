library(dplyr)
library(readr)

mutate_columns <- function(df) {
  df %>% mutate(
    dataset = snakemake@params[["dataset"]],
    ref = snakemake@params[["ref"]],
    region = snakemake@params[["region"]],
    seed = snakemake@params[["seed"]],
    method = snakemake@params[["method"]],
    printref = snakemake@params[["printref"]]
  )
}

reformat <- function(key) {
  read_tsv(snakemake@input[[key]]) %>% 
    mutate_columns %>% 
    write_tsv(snakemake@output[[key]])
}

reformat('sensspec')
reformat('bench')
