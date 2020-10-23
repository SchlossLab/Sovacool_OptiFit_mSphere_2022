library(dplyr)
readr::read_tsv(snakemake@input[["count"]]) %>%
  mutate(
    dummyRefGroup = total,
    Representative_Sequence = gsub("[\\.-]", "_", Representative_Sequence)
  ) %>%
  select(Representative_Sequence, dummyRefGroup) %>%
  readr::write_tsv(snakemake@output[["count"]])
