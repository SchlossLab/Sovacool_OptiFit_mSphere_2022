library(here)
library(tidyverse)

rel_diff <- function(final, init, percent = TRUE) {
  mult <- if (isTRUE(percent)) 100 else 1
  return((final - init) / init * mult)
}
coeff_var <- function(x) {
  return(sd(x) / mean(x))
}


dat <- read_tsv(here("results", "summarized.tsv"))
agg <- read_tsv(here("results", "aggregated.tsv"))
################################################################################
# de novo datasets
opticlust_mcc <- agg %>%
  filter(
    method == "de_novo",
    tool == "mothur"
  ) %>%
  pull(mcc) %>%
  median()
opticlust_sec <- agg %>%
  filter(
    method == "de_novo",
    tool == "mothur"
  ) %>%
  pull(sec) %>%
  median()
opticlust_mem <- agg %>%
  filter(
    strategy == "de_novo",
    tool == "mothur"
  ) %>%
  pull(mem_gb) %>%
  median()
dn_vsearch_mcc <- agg %>%
  filter(strategy == "de_novo", tool == "vsearch") %>%
  pull(mcc) %>%
  median()
dn_vsearch_sec <- agg %>%
  filter(strategy == "de_novo", tool == "vsearch") %>%
  pull(sec) %>%
  median()
mcc_opticlust_vs_vsearch <- rel_diff(opticlust_mcc, dn_vsearch_mcc)
sec_opticlust_vs_vsearch <- abs(rel_diff(dn_vsearch_sec, opticlust_sec))

################################################################################
# de novo ref dbs

dn_dbs <- read_tsv("subworkflows/2_fit_reference_db/results/denovo_dbs.tsv") %>%
  group_by(ref) %>%
  summarize(med_mcc = median(mcc)) %>%
  full_join(read_tsv(here(
    "subworkflows", "0_prep_db", "data",
    "seq_counts.tsv"
  )),
  by = "ref"
  ) %>%
  mutate(refname = case_when(
    ref == "gg" ~ "Greengenes",
    TRUE ~ toupper(ref)
  ))

################################################################################
# ref db open
open_fit_db_mcc <- agg %>%
  filter(
    method == "open",
    strategy == "database",
    tool == "mothur"
  ) %>%
  pull(mcc) %>%
  median()

mcc_open_fit_db_vs_clust <- rel_diff(opticlust_mcc, open_fit_db_mcc)

open_fit_gg_mcc <- agg %>%
  filter(
    method == "open",
    strategy == "database",
    tool == "mothur",
    ref == "gg"
  ) %>%
  pull(mcc) %>%
  median()

open_fit_silva_mcc <- agg %>%
  filter(
    method == "open",
    strategy == "database",
    tool == "mothur",
    ref == "silva"
  ) %>%
  pull(mcc) %>%
  median()

open_fit_rdp_mcc <- agg %>%
  filter(
    method == "open",
    strategy == "database",
    tool == "mothur",
    ref == "rdp"
  ) %>%
  pull(mcc) %>%
  median()

open_vsearch_mcc <- agg %>%
  filter(
    method == "open",
    strategy == "database",
    tool == "vsearch",
    ref == "gg"
  ) %>%
  pull(mcc) %>%
  median()
mcc_open_fit_db_vs_vsearch <- rel_diff(open_fit_gg_mcc, open_vsearch_mcc)

open_vsearch_sec <- agg %>%
  filter(
    method == "open",
    strategy == "database",
    tool == "vsearch",
    ref == "gg"
  ) %>%
  pull(sec) %>%
  median()
open_fit_db_sec <- agg %>%
  filter(
    method == "open",
    strategy == "database",
    tool == "mothur"
  ) %>%
  pull(sec) %>%
  median()
sec_vsearch_vs_open_fit_db <- rel_diff(open_vsearch_sec, open_fit_db_sec)
sec_opticlust_vs_open_fit_db <- rel_diff(opticlust_sec, open_fit_db_sec) %>% abs()

# human dataset to silva
open_fit_silva_human_sec <- agg %>%
  filter(
    method == "open",
    strategy == "database",
    tool == "mothur",
    ref == "silva",
    dataset == "human"
  ) %>%
  pull(sec) %>%
  median()
closed_fit_silva_human_sec <- agg %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur",
    ref == "silva",
    dataset == "human"
  ) %>%
  pull(sec) %>%
  median()
opticlust_human_sec <- agg %>%
  filter(
    method == "de_novo",
    tool == "mothur",
    dataset == "human"
  ) %>%
  pull(sec) %>%
  median()

################################################################################
# ref db closed
closed_fit_db_mcc <- agg %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur"
  ) %>%
  pull(mcc) %>%
  median()
mcc_closed_fit_db_vs_clust <- rel_diff(closed_fit_db_mcc, opticlust_mcc) %>% abs()

closed_fit_gg_mcc <- agg %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur",
    ref == "gg"
  ) %>%
  pull(mcc) %>%
  median()
closed_fit_silva_mcc <- agg %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur",
    ref == "silva"
  ) %>%
  pull(mcc) %>%
  median()
closed_fit_rdp_mcc <- agg %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur",
    ref == "rdp"
  ) %>%
  pull(mcc) %>%
  median()

frac_fit_db <- agg %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur"
  ) %>%
  pull(fraction_mapped) %>%
  median() * 100
frac_fit_gg <- agg %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur",
    ref == "gg"
  ) %>%
  pull(fraction_mapped) %>%
  median() * 100
frac_fit_silva <- agg %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur",
    ref == "silva"
  ) %>%
  pull(fraction_mapped) %>%
  median() * 100
frac_fit_rdp <- agg %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur",
    ref == "rdp"
  ) %>%
  pull(fraction_mapped) %>%
  median() * 100
frac_vsearch <- agg %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "vsearch",
    ref == "gg"
  ) %>%
  pull(fraction_mapped) %>%
  median() * 100
frac_vsearch_vs_fit <- rel_diff(frac_vsearch, frac_fit_gg)

closed_fit_db_sec <- agg %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur"
  ) %>%
  pull(sec) %>%
  median()
sec_closed_fit_db_vs_clust <- rel_diff(closed_fit_db_sec, opticlust_sec) %>% abs()

closed_vsearch_sec <- agg %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "vsearch"
  ) %>%
  pull(sec) %>%
  median()
sec_closed_fit_db_vs_vsearch <- rel_diff(closed_fit_db_sec, closed_vsearch_sec) %>% abs()

closed_vsearch_mcc <- agg %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "vsearch"
  ) %>%
  pull(mcc) %>%
  median()
mcc_closed_fit_db_vs_vsearch <- rel_diff(closed_fit_db_mcc, closed_vsearch_mcc) %>% abs()

################################################################################
# fit split
cv_fit_split_mcc <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "simple",
    ref_frac == 0.5
  ) %>%
  pull(mcc) %>%
  coeff_var()

frac_fit_split <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "simple",
    method == "closed",
    ref_frac == 0.5
  ) %>%
  pull(fraction_mapped) %>%
  median() * 100

closed_fit_split_sec <- agg %>%
  filter(
    method == "closed",
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "simple",
    ref_frac == 0.5
  ) %>%
  pull(sec) %>%
  median()
open_fit_split_sec <- agg %>%
  filter(
    method == "open",
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "simple",
    ref_frac == 0.5
  ) %>%
  pull(sec) %>%
  median()

sec_closed_fit_split_vs_clust <- rel_diff(opticlust_sec, closed_fit_split_sec) %>% abs()
sec_open_fit_split_vs_clust <- rel_diff(opticlust_sec, open_fit_split_sec) %>% abs()
sec_open_fit_split_vs_db <- rel_diff(open_fit_db_sec, open_fit_split_sec) %>% abs()
sec_closed_fit_split_vs_db <- rel_diff(closed_fit_db_sec, closed_fit_split_sec) %>% abs()

closed_fit_split_mem <- agg %>%
  filter(
    method == "closed",
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "simple",
    ref_frac == 0.5
  ) %>%
  pull(mem_gb) %>%
  median()
open_fit_split_mem <- agg %>%
  filter(
    method == "open",
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "simple",
    ref_frac == 0.5
  ) %>%
  pull(mem_gb) %>%
  median()
mem_closed_fit_split_vs_clust <- rel_diff(closed_fit_split_mem, opticlust_mem)
mem_open_fit_split_vs_clust <- rel_diff(open_fit_split_mem, opticlust_mem)

cv_fit_split_mcc_human_simple <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    dataset == "human",
    ref_weight == "simple"
  ) %>%
  pull(mcc) %>%
  coeff_var()

cv_fit_split_mem_human_simple <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    dataset == "human",
    ref_weight == "simple"
  ) %>%
  pull(mem_gb) %>%
  coeff_var()

sec_fit_split_human_simple_1 <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    dataset == "human",
    ref_weight == "simple",
    ref_frac == 0.1
  ) %>%
  pull(sec) %>%
  median()

sec_fit_split_human_simple_9 <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    dataset == "human",
    ref_weight == "simple",
    ref_frac == 0.9
  ) %>%
  pull(sec) %>%
  median()

frac_fit_split_human_simple_1 <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    dataset == "human",
    ref_weight == "simple",
    method == "closed",
    ref_frac == 0.1
  ) %>%
  pull(fraction_mapped) %>%
  median()

frac_fit_split_human_simple_9 <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    dataset == "human",
    ref_weight == "simple",
    method == "closed",
    ref_frac == 0.9
  ) %>%
  pull(fraction_mapped) %>%
  median()

mcc_fit_split_simple <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "simple",
    ref_frac == 0.5
  ) %>%
  pull(mcc) %>%
  median()
mcc_opticlust_vs_fit_split_simple <- rel_diff(opticlust_mcc, mcc_fit_split_simple)

#####
# fit split at ref_frac 0.5

mcc_fit_split_abun <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "abundance",
    ref_frac == 0.5
  ) %>%
  pull(mcc) %>%
  median()

mcc_fit_split_dist <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "distance",
    ref_frac == 0.5
  ) %>%
  pull(mcc) %>%
  median()

frac_fit_split_simple <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "simple",
    method == "closed",
    ref_frac == 0.5
  ) %>%
  pull(fraction_mapped) %>%
  median()

frac_fit_split_abun <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "abundance",
    method == "closed",
    ref_frac == 0.5
  ) %>%
  pull(fraction_mapped) %>%
  median()

frac_fit_split_dist <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "distance",
    method == "closed",
    ref_frac == 0.5
  ) %>%
  pull(fraction_mapped) %>%
  median()

sec_fit_split_simple <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "simple",
    ref_frac == 0.5
  ) %>%
  pull(sec) %>%
  median()

sec_fit_split_abun <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "abundance",
    ref_frac == 0.5
  ) %>%
  pull(sec) %>%
  median()

sec_fit_split_dist <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "distance",
    ref_frac == 0.5
  ) %>%
  pull(sec) %>%
  median()

##########

frac_fit_split_1 <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    method == "closed",
    ref_frac == 1
  ) %>%
  pull(fraction_mapped) %>%
  median() * 100

frac_fit_split_9 <- agg %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    method == "closed",
    ref_frac == 0.9
  ) %>%
  pull(fraction_mapped) %>%
  median() * 100


################################################################################
# save results
save.image(file = here("results", "stats.RData"))
