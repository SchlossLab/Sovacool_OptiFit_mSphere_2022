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
