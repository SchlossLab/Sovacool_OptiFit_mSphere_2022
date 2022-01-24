library(tidyverse)
library(here)
dat_agg <- read_tsv(here('results', 'aggregated.tsv'))
dat_agg %>% 
  summarise(secs = sum(sec)) %>% 
  mutate(mins = secs/60, hrs = mins/60, days = hrs/24)
# A tibble: 1 × 4
#       secs   mins   hrs  days
#      <dbl>  <dbl> <dbl> <dbl>
# 1 2130410. 35507.  592.  24.7
dat <- data.frame(dx = c('normal', 'cancer', 'normal', 'cancer'),
                  condition = c('A', 'A', 'B', 'B'),
                  otu1 = c(1, 0, 0, 1),
                  otu2 = c(0, 1, 0, 1)
                  )
conditions <- dat %>% pull('condition')
dat <- dat %>% select(-conditions)
run_ml(dat, 'rf', outcome_colname = 'dx', groups = conditions, ...)