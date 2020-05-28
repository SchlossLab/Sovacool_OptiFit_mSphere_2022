plots.Rmd
================
Kelly Sovacool
5/27/2020

``` r
library(here)
library(tidyverse)
theme_set(theme_classic())
```

# Fitting datasets to a reference database vs *de novo* clustering

``` r
sensspec <- read_tsv(here('subworkflows/2_fit_reference_db/results/merged.sensspec')) %>% 
  mutate(method=fct_relevel(method, c("open", "closed", "de_novo")))
sensspec %>% group_by(dataset, method) %>% 
  ggplot(aes(x=method, y=mcc, fill=dataset)) +
  geom_jitter()
```

![](figures/fit_db_sensspec-1.png)<!-- -->

``` r
benchmarks <- read_tsv(here('subworkflows/2_fit_reference_db/results/benchmarks.txt')) %>% 
  mutate(method=fct_relevel(method, c("open", "closed", "de_novo")))
benchmarks %>% group_by(dataset, method) %>% 
  ggplot(aes(x=method, y=s, fill=dataset)) +
  geom_jitter() +
  labs(y='seconds')
```

![](figures/fit_db_benchmarks-1.png)<!-- -->
