library(tidyverse)
library(glue)
library(here)

mutate_perf <- function(dat) {
  dat %>%
    mutate(
      mem_mb = max_rss,
      mem_gb = mem_mb / 1024
    ) %>%
    rename(
      sec = s
    )
}
#
# opticlust <- read_tsv(snakemake@input[['opticlust']])
# optifit_db <- read_tsv(snakemake@input[['optifit_db']])
# optifit_split <- read_tsv(snakemake@input[['optifit_split']])
# vsearch <- read_tsv(snakemake@input[['vsearch']])

opticlust <- read_tsv(here("subworkflows/1_prep_samples/results/opticlust_results.tsv")) %>%
  mutate_perf() %>%
  mutate(strategy = method)
optifit_db <- read_tsv(here("subworkflows/2_fit_reference_db/results/optifit_dbs_results.tsv")) %>%
  mutate_perf() %>%
  mutate(strategy = "database")
optifit_split <- read_tsv(here("subworkflows/3_fit_sample_split/results/optifit_split_results.tsv")) %>%
  mutate_perf() %>%
  mutate(strategy = "self-split")
vsearch <- read_tsv(here("subworkflows/4_vsearch/results/vsearch_results.tsv")) %>%
  mutate_perf() %>%
  mutate(strategy = case_when(
    method == "de_novo" ~ method,
    TRUE ~ "database"
  ))
results_agg <- list(opticlust, optifit_db, optifit_split, vsearch) %>%
  reduce(full_join)

results_sum <- results_agg %>%
  group_by(tool, strategy, method, dataset, ref, ref_frac, ref_weight) %>%
  summarize(
    n = n(),
    mcc_median = median(mcc), # TODO: tidy way to avoid this repetitiveness?
    sec_median = median(sec),
    mem_gb_median = median(mem_gb),
    frac_map_median = median(fraction_mapped)
  )
#
# write_tsv(results_agg, snakemake@output[['agg']])
# write_tsv(results_sum, snakemake@output[['sum']])
write_tsv(results_agg, "results/aggregated.tsv")
write_tsv(results_sum, "results/summarized.tsv")
