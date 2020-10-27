library(dplyr)
readr::read_tsv(snakemake@input[["count"]]) %>%
  mutate(
    dummyRefGroup = total, # vsearch doesn't support dots or hyphens in seq names
    Representative_Sequence = gsub("[\\.-]", "_", Representative_Sequence)
  ) %>%
  select(Representative_Sequence, total, dummyRefGroup) %>%
  readr::write_tsv(snakemake@output[["count"]])
