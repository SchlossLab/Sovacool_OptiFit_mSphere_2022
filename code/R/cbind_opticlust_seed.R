source(snakemake@input[["fcns"]])
log_smk()

mutate_columns <- function(df_orig) {
  df_new <- df_orig %>% mutate(
    dataset = snakemake@wildcards[["dataset"]],
    ref = NA,
    region = NA,
    seed = snakemake@wildcards[["seed"]],
    method = "de_novo",
    printref = NA,
    tool = "mothur"
  )
  if (all(c("ref_weight", "ref_frac") %in% names(snakemake@wildcards))) {
    df_new <- df_new %>% mutate(
      ref_weight = snakemake@wildcards[["ref_weight"]],
      ref_frac = snakemake@wildcards[["ref_frac"]],
      sample_frac = NA
    )
  }
  return(df_new)
}

cbind_all("tsv")
