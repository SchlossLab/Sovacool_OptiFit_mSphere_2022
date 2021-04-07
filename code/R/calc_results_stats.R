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

################################################################################
# de novo
opticlust_mcc <- dat %>%
  filter(
    method == "de_novo",
    tool == "mothur"
  ) %>%
  pull(mcc_median) %>%
  median()
opticlust_sec <- dat %>%
  filter(
    method == "de_novo",
    tool == "mothur"
  ) %>%
  pull(sec_median) %>%
  median()
dn_vsearch_mcc <-
  dat %>%
  filter(strategy == "de_novo", tool == "vsearch") %>%
  pull(mcc_median) %>%
  median()
dn_vsearch_sec <-
  dat %>%
  filter(strategy == "de_novo", tool == "vsearch") %>%
  pull(sec_median) %>%
  median()
mcc_opticlust_vs_vsearch <- rel_diff(opticlust_mcc, dn_vsearch_mcc)
sec_opticlust_vs_vsearch <- abs(rel_diff(dn_vsearch_sec, opticlust_sec))

################################################################################
# ref db open
open_fit_db_mcc <- dat %>%
  filter(
    method == "open",
    strategy == "database",
    tool == "mothur"
  ) %>%
  pull(mcc_median) %>%
  median()

mcc_open_fit_db_vs_clust <- rel_diff(opticlust_mcc, open_fit_db_mcc)

open_fit_gg_mcc <- dat %>%
  filter(
    method == "open",
    strategy == "database",
    tool == "mothur",
    ref == "gg"
  ) %>%
  pull(mcc_median) %>%
  median()

open_fit_silva_mcc <- dat %>%
  filter(
    method == "open",
    strategy == "database",
    tool == "mothur",
    ref == "silva"
  ) %>%
  pull(mcc_median) %>%
  median()

open_fit_rdp_mcc <- dat %>%
  filter(
    method == "open",
    strategy == "database",
    tool == "mothur",
    ref == "rdp"
  ) %>%
  pull(mcc_median) %>%
  median()

open_vsearch_mcc <- dat %>%
  filter(
    method == "open",
    strategy == "database",
    tool == "vsearch",
    ref == "gg"
  ) %>%
  pull(mcc_median) %>%
  median()
mcc_open_fit_db_vs_vsearch <- rel_diff(open_fit_gg_mcc, open_vsearch_mcc)

open_vsearch_sec <- dat %>%
  filter(
    method == "open",
    strategy == "database",
    tool == "vsearch",
    ref == "gg"
  ) %>%
  pull(sec_median) %>%
  median()
open_fit_db_sec <- dat %>%
  filter(
    method == "open",
    strategy == "database",
    tool == "mothur"
  ) %>%
  pull(sec_median) %>%
  median()
sec_vsearch_vs_open_fit_db <- rel_diff(open_vsearch_sec, open_fit_db_sec)
sec_opticlust_vs_open_fit_db <- rel_diff(opticlust_sec, open_fit_db_sec) %>% abs()

################################################################################
# ref db closed
closed_fit_db_mcc <- dat %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur"
  ) %>%
  pull(mcc_median) %>%
  median()
mcc_closed_fit_db_vs_clust <- rel_diff(closed_fit_db_mcc, opticlust_mcc) %>% abs()

closed_fit_gg_mcc <- dat %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur",
    ref == "gg"
  ) %>%
  pull(mcc_median) %>%
  median()
closed_fit_silva_mcc <- dat %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur",
    ref == "silva"
  ) %>%
  pull(mcc_median) %>%
  median()
closed_fit_rdp_mcc <- dat %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur",
    ref == "rdp"
  ) %>%
  pull(mcc_median) %>%
  median()

frac_fit_db <- dat %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur"
  ) %>%
  pull(frac_map_median) %>%
  median() * 100
frac_fit_gg <- dat %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur",
    ref == "gg"
  ) %>%
  pull(frac_map_median) %>%
  median() * 100
frac_fit_silva <- dat %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur",
    ref == "silva"
  ) %>%
  pull(frac_map_median) %>%
  median() * 100
frac_fit_rdp <- dat %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur",
    ref == "rdp"
  ) %>%
  pull(frac_map_median) %>%
  median() * 100
frac_vsearch <- dat %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "vsearch",
    ref == "gg"
  ) %>%
  pull(frac_map_median) %>%
  median() * 100
frac_vsearch_vs_fit <- rel_diff(frac_vsearch, frac_fit_gg)

closed_fit_db_sec <- dat %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "mothur"
  ) %>%
  pull(sec_median) %>%
  median()
sec_closed_fit_db_vs_clust <- rel_diff(closed_fit_db_sec, opticlust_sec) %>% abs()
closed_vsearch_sec <- dat %>%
  filter(
    method == "closed",
    strategy == "database",
    tool == "vsearch"
  ) %>%
  pull(sec_median) %>%
  median()
sec_closed_fit_db_vs_vsearch <- rel_diff(closed_fit_db_sec, closed_vsearch_sec) %>% abs()


################################################################################
# fit split
cv_fit_split_mcc <- coeff_var(dat %>% filter(
  strategy == "self-split",
  tool == "mothur",
  ref_weight == 'simple',
  ref_frac == 0.5
) %>%
  pull(mcc_median))
cv_fit_split_mcc_human <- coeff_var(dat %>% filter(
  strategy == "self-split",
  tool == "mothur",
  dataset == "human",
  ref_weight == 'simple',
  ref_frac == 0.5
) %>%
  pull(mcc_median))
closed_fit_split_sec <- dat %>%
  filter(
    method == "closed",
    strategy == "self-split",
    tool == "mothur",
    ref_weight == 'simple',
    ref_frac == 0.5
  ) %>%
  pull(sec_median) %>%
  median()
sec_closed_fit_split_vs_clust <- rel_diff(closed_fit_split_sec, opticlust_sec) %>% abs()
open_fit_split_sec <- dat %>%
  filter(
    method == "open",
    strategy == "self-split",
    tool == "mothur",
    ref_weight == 'simple',
    ref_frac == 0.5
  ) %>%
  pull(sec_median) %>%
  median()
sec_open_fit_split_vs_clust <- rel_diff(open_fit_split_sec, opticlust_sec) %>% abs()
sec_open_fit_split_vs_db <- rel_diff(open_fit_split_sec, open_fit_db_sec) %>% abs()
sec_closed_fit_split_vs_db <- rel_diff(closed_fit_split_sec, closed_fit_db_sec) %>% abs()

mcc_fit_split_simple <- dat %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "simple",
    ref_frac == 0.5
  ) %>%
  pull(mcc_median) %>%
  median()
mcc_opticlust_vs_fit_split_simple <- rel_diff(opticlust_mcc, mcc_fit_split_simple)

mcc_fit_split_abun <- dat %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "abundance",
    ref_frac == 0.5
  ) %>%
  pull(mcc_median) %>%
  median()

mcc_fit_split_dist <- dat %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    ref_weight == "distance",
    ref_frac == 0.5
  ) %>%
  pull(mcc_median) %>%
  median()

frac_fit_split_0.1 <- dat %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    method == "closed",
    ref_frac == 0.1
  ) %>%
  pull(frac_map_median) %>%
  median() * 100


frac_fit_split_0.8 <- dat %>%
  filter(
    strategy == "self-split",
    tool == "mothur",
    method == "closed",
    ref_frac == 0.8
  ) %>%
  pull(frac_map_median) %>%
  median() * 100


################################################################################
# save results
save.image(file = here("results", "stats.RData"))
