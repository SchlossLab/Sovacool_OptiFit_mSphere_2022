library(dplyr)
library(purrr)
library(readr)

#' Write all messages to a log file if specified in the snakemake workflow
#'
#' @export
#'
log_smk <- function() {
  if (!is.null(snakemake@log)) {
    log_filepath <- snakemake@log[1][[1]]
    log <- file(log_filepath, open = "wt")
    sink(log, append = TRUE)
    sink(log, append = TRUE, type = "message")
  }
}

#' Merge rows of all tsv files to one file
#'
#' @param key name of list of tsv files in snakemake input field
#'
#' @return (None) writes tsv file
#' @export
#'
merge_results <- function(key) {
  snakemake@input[[key]] %>%
    map(read_tsv) %>%
    bind_rows() %>%
    write_tsv(snakemake@output[[key]])
}

mutate_columns <- function(df_orig) {
  df_new <- df_orig %>% mutate(
    dataset = snakemake@params[["dataset"]],
    ref = snakemake@params[["ref"]],
    region = snakemake@params[["region"]],
    seed = snakemake@params[["seed"]],
    method = snakemake@params[["method"]],
    printref = snakemake@params[["printref"]]
  )
  if (all(c("ref_weight", "ref_frac", "sample_frac") %in% names(snakemake@params))) {
    df_new <- df_new %>% mutate(
      ref_weight = snakemake@params[["ref_weight"]],
      ref_frac = snakemake@params[["ref_frac"]],
      sample_frac = snakemake@params[["sample_frac"]]
    )
  }
  return(df_new)
}

reformat <- function(key) {
  read_tsv(snakemake@input[[key]]) %>%
    mutate_columns() %>%
    write_tsv(snakemake@output[[key]])
}