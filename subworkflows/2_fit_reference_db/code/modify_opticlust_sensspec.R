library(dplyr)
library(readr)
sensspec <- read_tsv(snakemake@input[["txt"]]) %>%
  mutate(
    dataset = snakemake@params[["dataset"]],
    ref = snakemake@params[["ref"]],
    region = snakemake@params[["region"]],
    seed = snakemake@params[["seed"]],
    method = "de_novo",
    printref = NA,
    iter = NA,
    numotus = NA
  )
write_tsv(sensspec, snakemake@output[["txt"]])
