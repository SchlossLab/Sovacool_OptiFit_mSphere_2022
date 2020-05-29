Exploratory Plots
================
5/27/2020

``` r
library(here)
library(tidyverse)
theme_set(theme_classic())
```

# Fitting datasets to a reference database vs *de novo* clustering

``` r
relevel_method <- function(df) {
  df %>%
    mutate(method = fct_relevel(method, c("open", "closed", "de_novo")))
}
plot_jitter <- function(df, y) {
  df %>% group_by(dataset, method) %>% 
      ggplot(aes(x=method, y={{ y }}, color=dataset)) +
  geom_jitter()
}
```

``` r
sensspec <- read_tsv(here('subworkflows/2_fit_reference_db/results/sensspec.txt')) %>% 
  relevel_method()
sensspec %>% 
  plot_jitter(mcc)
```

![](figures/fit_db_sensspec-1.png)<!-- -->

``` r
benchmarks <- read_tsv(here('subworkflows/2_fit_reference_db/results/benchmarks.txt')) %>% 
  relevel_method
benchmarks %>% 
  plot_jitter(s) +
  labs(y='seconds')
```

![](figures/fit_db_benchmarks-1.png)<!-- -->

``` r
fractions <- read_tsv(here('subworkflows/2_fit_reference_db/results/fraction_reads_mapped.txt'))
fractions %>% 
  ggplot(aes(x=dataset, y=fraction_mapped, color=dataset)) +
  geom_jitter() +
  ylim(0, 1)
```

![](figures/fraction_reads_mapped-1.png)<!-- -->
