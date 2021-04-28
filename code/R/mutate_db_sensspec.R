source(snakemake@input[["fcns"]])
log_smk()

na_if_not <- function(key) {
  ifelse(key %in% names(snakemake@wildcards),
    snakemake@wildcards[[key]],
    NA
  )
}

mutate_columns <- function(df_orig) {
  df_orig %>% mutate(
    ref = na_if_not("ref"),
    dataset_filter = na_if_not("dataset"),
    region = na_if_not("region"),
    seed = na_if_not("seed"),
    method = "de_novo",
    tool = "mothur"
  )
}

read_tsv(snakemake@input[['tsv']]) %>%
    mutate_colums() %>%
    write_tsv(snakemake@output[['tsv']])