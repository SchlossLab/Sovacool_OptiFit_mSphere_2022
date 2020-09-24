source(snakemake@input[["fcns"]])

log_smk()
snakemake@input[["tsv"]] %>%
  purrr::map(readr::read_tsv) %>%
  purrr::reduce(dplyr::full_join) %>%
  readr::write_tsv(snakemake@output[["tsv"]])
