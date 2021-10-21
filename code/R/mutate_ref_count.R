library(dplyr)
readr::read_tsv(snakemake@input[["count_table"]]) %>%
  mutate(dummyRefGroup = total) %>%
  select(Representative_Sequence, total, dummyRefGroup) %>%
  readr::write_tsv(snakemake@output[["count_table"]])
