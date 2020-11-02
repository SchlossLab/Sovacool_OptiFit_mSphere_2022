source(snakemake@input[["fcns"]])
log_smk()

na_if_not <- function(key) {
    ifelse(key %in% names(snakemake@wildcards),
           snakemake@wildcards[[key]],
           NA)
}

mutate_columns <- function(df_orig) {
  df_orig %>% mutate(
    dataset = na_if_not("dataset"),
    ref = na_if_not("ref"),
    region = na_if_not("region"),
    seed = na_if_not("seed"),
    method = na_if_not("method"),
    printref = na_if_not("printref"),
    ref_weight = na_if_not("ref_weight"),
    ref_frac = na_if_not("ref_frac"),
    sample_frac = na_if_not("sample_frac"),
    tool = "mothur"
    )
}

rbind_all('tsv')