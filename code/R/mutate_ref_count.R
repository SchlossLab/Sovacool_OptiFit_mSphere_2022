readr::read_tsv(snakemake@input[["count"]]) %>%
  mutate(dummyRefGroup = total) %>%
  readr::write_tsv(snakemake@output[["count"]])
