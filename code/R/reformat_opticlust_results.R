source(snakemake@input[["fcns"]])
log_smk()

mutate_columns <- function(df_orig) {
  df_new <- df_orig %>% mutate(
    dataset = snakemake@params[["dataset"]],
    ref = NA,
    region = NA,
    seed = snakemake@params[["seed"]],
    method = "de_novo",
    printref = NA,
    tool = "mothur"
  )
  if (all(c("ref_weight", "ref_frac") %in% names(snakemake@params))) {
    df_new <- df_new %>% mutate(
      ref_weight = snakemake@params[["ref_weight"]],
      ref_frac = snakemake@params[["ref_frac"]],
      sample_frac = NA
    )
  }
  return(df_new)
}

reformat <- function(key) {
  read_tsv(snakemake@input[[key]]) %>%
    mutate_columns() %>%
    write_tsv(snakemake@output[[key]])
}

reformat("sensspec")
reformat("bench")
reformat("div")
