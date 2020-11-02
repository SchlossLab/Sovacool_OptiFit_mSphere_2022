source(snakemake@input[["fcns"]])
log_smk()

mutate_columns <- function(df_orig) {
  df_new <- df_orig %>% mutate(
    dataset = snakemake@wildcards[["dataset"]],
    ref = snakemake@wildcards[["ref"]],
    region = snakemake@wildcards[["region"]],
    seed = snakemake@wildcards[["seed"]],
    method = snakemake@wildcards[["method"]],
    printref = snakemake@wildcards[["printref"]],
    tool = "mothur"
  )
  if (all(c("ref_weight", "ref_frac", "sample_frac") %in% names(snakemake@wildcards))) {
    df_new <- df_new %>% mutate(
      ref_weight = snakemake@wildcards[["ref_weight"]],
      ref_frac = snakemake@wildcards[["ref_frac"]],
      sample_frac = snakemake@wildcards[["sample_frac"]]
    )
  }
  return(df_new)
}

rbind_all('tsv')