library(here)
library(tidyverse)

rel_diff <- function(final, init, percent = TRUE) {
  mult <- if (isTRUE(percent)) 100 else 1
  return((final - init) / init * mult)
}


dat <- read_tsv(here('results', 'summarized.tsv'))

# de novo
opticlust_mcc <- dat %>% filter(method == 'de_novo', 
                                tool == 'mothur') %>% 
  pull(mcc_median) %>% median()
opticlust_sec <- dat %>% filter(method == 'de_novo', 
                                tool == 'mothur') %>% 
  pull(sec_median) %>% median()
dn_vsearch_mcc <-
  dat %>% filter(strategy == 'de_novo', tool == 'vsearch') %>% 
  pull(mcc_median) %>% median()
dn_vsearch_sec <-
  dat %>% filter(strategy == 'de_novo', tool == 'vsearch') %>% 
  pull(sec_median) %>% median()
mcc_opticlust_vs_vsearch <- rel_diff(opticlust_mcc, dn_vsearch_mcc)
sec_opticlust_vs_vsearch <- abs(rel_diff(dn_vsearch_sec, opticlust_sec))

# ref db open 
open_fit_mcc <- dat %>% filter(method == 'open', 
                               strategy == 'database', 
                               tool == 'mothur') %>% 
  pull(mcc_median) %>% median()
mcc_open_fit_db_vs_clust <- rel_diff(opticlust_mcc, open_fit_mcc)

open_fit_gg_mcc <- dat %>% filter(method == 'open',
                                  strategy == 'database',
                                  tool == 'mothur',
                                  ref == 'gg') %>% pull(mcc_median) %>% median()
open_vsearch_mcc <- dat %>% filter(method == 'open',
                                   strategy == 'database',
                                   tool == 'vsearch',
                                   ref == 'gg') %>% pull(mcc_median) %>% median()
mcc_open_fit_db_vs_vsearch <- rel_diff(open_fit_gg_mcc, open_vsearch_mcc)

open_vsearch_sec <- dat %>% filter(method == 'open',
                                   strategy == 'database',
                                   tool == 'vsearch',
                                   ref == 'gg') %>% pull(sec_median) %>% median()
open_fit_sec <- dat %>% filter(method == 'open', 
                               strategy == 'database', 
                               tool == 'mothur') %>%
  pull(sec_median) %>% median()
sec_vsearch_vs_open_fit_db <- rel_diff(open_vsearch_sec, open_fit_sec)
sec_opticlust_vs_open_fit_db <- rel_diff(open_fit_sec, opticlust_sec)

# ref db closed
closed_fit_db_mcc <- dat %>% filter(method == 'closed',
                                    strategy == 'database',
                                    tool == 'mothur') %>% 
  pull(mcc_median) %>% median()
mcc_closed_fit_db_vs_clust <- rel_diff(closed_fit_db_mcc, opticlust_mcc) %>% abs()

closed_fit_gg_mcc <- dat %>% filter(method == 'closed',
                                    strategy == 'database',
                                    tool == 'mothur',
                                    ref == 'gg') %>% 
  pull(mcc_median) %>% median()
closed_fit_silva_mcc <- dat %>% filter(method == 'closed',
                                       strategy == 'database',
                                       tool == 'mothur',
                                       ref == 'silva') %>% 
  pull(mcc_median) %>% median()
closed_fit_rdp_mcc <- dat %>% filter(method == 'closed',
                                     strategy == 'database',
                                     tool == 'mothur',
                                     ref == 'rdp') %>% 
  pull(mcc_median) %>% median()

frac_fit_db <-  dat %>% filter(method == 'closed',
                               strategy == 'database',
                               tool == 'mothur') %>%
  pull(frac_map_median) %>% median()
frac_fit_gg <-  dat %>% filter(method == 'closed',
                               strategy == 'database',
                               tool == 'mothur',
                               ref == 'gg') %>% 
  pull(frac_map_median) %>% median()
frac_fit_silva <-  dat %>% filter(method == 'closed',
                                  strategy == 'database',
                                  tool == 'mothur',
                                  ref == 'silva') %>% 
  pull(frac_map_median) %>% median()
frac_fit_rdp <-  dat %>% filter(method == 'closed',
                                strategy == 'database',
                                tool == 'mothur',
                                ref == 'rdp') %>% 
  pull(frac_map_median) %>% median()
frac_vsearch <-  dat %>% filter(method == 'closed',
                                strategy == 'database',
                                tool == 'vsearch',
                                ref == 'gg') %>% 
  pull(frac_map_median) %>% median()
frac_vsearch_vs_fit <- rel_diff(frac_vsearch, frac_fit_gg)

closed_fit_db_sec <- dat %>% filter(method == 'closed',
                                    strategy == 'database',
                                    tool == 'mothur') %>% 
  pull(sec_median) %>% median()
sec_closed_fit_db_vs_clust <- rel_diff(closed_fit_db_sec, opticlust_sec)
closed_vsearch_sec <- dat %>% filter(method == 'closed',
                                     strategy == 'database',
                                     tool == 'vsearch') %>% 
  pull(sec_median) %>% median()
sec_closed_fit_db_vs_vsearch <- rel_diff(closed_fit_db_sec, closed_vsearch_sec)

# save results
save.image(file = here('results', 'stats.RData'))