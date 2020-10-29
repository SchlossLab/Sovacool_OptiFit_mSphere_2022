library(dplyr)
library(purrr)
library(readr)

#' Write all messages to a log file if specified in the snakemake workflow
#'
#' @export
#'
log_smk <- function() {
  if (exists('snakemake') & length(snakemake@log) != 0) {
    log <- file(snakemake@log[1][[1]], open = "wt")
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
rbind_all <- function(key) {
  snakemake@input[[key]] %>%
    map(read_tsv) %>%
    bind_rows() %>%
    write_tsv(snakemake@output[[key]])
}
