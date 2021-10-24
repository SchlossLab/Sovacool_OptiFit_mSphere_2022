library(dplyr)

if (exists("snakemake")) {
  input_filename <- snakemake@input[["count_table"]]
  output_filename <- snakemake@output[["count_table"]]
} else {
  stop("This script is only intended to be run from a Snakemake workflow")
  # input_filename <- "data/soil_silva_bact_v4/preclust_db/silva.bact_v4.filter.unique.precluster.pick.count_table"
  # output_filename <- "data/soil_silva_bact_v4/preclust_db/silva.bact_v4.filter.unique.precluster.pick.mod.count_table"
}


readr::read_tsv(input_filename) %>%
  mutate(dummyRefGroup = total) %>%
  select(Representative_Sequence, total, dummyRefGroup) %>%
  readr::write_tsv(output_filename)
