library(dplyr)
library(readr)

mutate_columns <- function(df_orig) {
  df_new <- df_orig %>% mutate(
    dataset = snakemake@params[["dataset"]],
    ref = snakemake@params[["ref"]],
    region = snakemake@params[["region"]],
    seed = snakemake@params[["seed"]],
    method = snakemake@params[["method"]],
    printref = snakemake@params[["printref"]]
  )
  if (all(c('ref_weight', 'ref_frac', 'sample_frac') %in% names(snakemake@params))) {
    df_new <- df_new %>% mutate(
      ref_weight = snakemake@params[['ref_weight']],
      ref_frac = snakemake@params[['ref_frac']],
      sample_frac = snakemake@params[['sample_frac']]
    )
  }
  return(df_new)
}

reformat <- function(key) {
  read_tsv(snakemake@input[[key]]) %>%
    mutate_columns() %>%
    write_tsv(snakemake@output[[key]])
}

log <- file(snakemake@log[[1]])
sink(log)
sink(log, append = TRUE, type = "message")
reformat("sensspec")
reformat("bench")
