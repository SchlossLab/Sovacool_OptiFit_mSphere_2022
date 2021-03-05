library(tidyverse)
library(glue)

mutate_perf <- function(dat) {
  dat %>%
    mutate(
      mem_mb = max_rss,
      mem_gb = mem_mb / 1024
    ) %>%
    rename(
      sec = s,
      num_otus = sobs
    )
}

opticlust <- read_tsv(snakemake@input[["opticlust"]]) %>%
  mutate_perf() %>%
  mutate(strategy = method)
optifit_db <- read_tsv(snakemake@input[["optifit_db"]]) %>%
  mutate_perf() %>%
  mutate(strategy = glue("database_{ref}"))
optifit_split <- read_tsv(snakemake@input[["optifit_split"]]) %>%
  mutate_perf() %>%
  mutate(strategy = "self-split")
vsearch <- read_tsv(snakemake@input[["vsearch"]]) %>%
  mutate_perf() %>%
  mutate(strategy = case_when(
    method == "de_novo" ~ method,
    TRUE ~ as.character(glue("database_{ref}"))
  ))

results_agg <- list(opticlust, optifit_db, optifit_split, vsearch) %>%
  reduce(full_join)

results_sum <- results_agg %>%
  group_by(tool, strategy, method, dataset, database, ref_frac, ref_weight) %>%
  summarize(
    n = n(),
    mcc_median = median(mcc), # TODO: tidy way to avoid this repetitiveness?
    sec_median = median(sec),
    mem_gb_median = median(mem_gb),
    frac_map_median = median(fraction_mapped)
  )

write_tsv(results_agg, snakemake@output[["agg"]])
write_tsv(results_sum, snakemake@output[["sum"]])
