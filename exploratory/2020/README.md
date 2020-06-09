Exploratory Plots
================
5/27/2020

``` r
library(here)
library(tidyverse)
theme_set(theme_classic())
color_palette <- RColorBrewer::brewer.pal(4, "Dark2")
dataset_colors <- c(
  human = color_palette[[3]],
  marine = color_palette[[1]],
  mouse = color_palette[[4]],
  soil = color_palette[[2]]
)
```

# Fitting datasets to a reference database vs *de novo* clustering

``` r
relevel_method <- function(df) {
  df %>%
    mutate(method = fct_relevel(method, c("open", "closed", "de_novo")))
}
```

``` r
sensspec <-
  read_tsv(here('subworkflows/2_fit_reference_db/results/sensspec.tsv')) %>%
  relevel_method()
sensspec %>%
  group_by(dataset, method) %>%
  ggplot(aes(x = method, y = mcc, color = dataset)) +
  geom_boxplot(alpha = 0.5) +
  scale_color_manual(values = dataset_colors) +
  #ylim(0, 1) +
  facet_grid(dataset ~ ref)
```

![](figures/fit_db_sensspec-1.png)<!-- -->

``` r
benchmarks <-
  read_tsv(here('subworkflows/2_fit_reference_db/results/benchmarks.tsv')) %>%
  relevel_method
benchmarks %>%
  group_by(dataset, method) %>%
  ggplot(aes(x = method, y = s, color = dataset)) +
  geom_boxplot() +
  scale_color_brewer(palette = "Dark2") +
  #scale_color_manual(values = dataset_colors) +  # bug in geom_boxplot with manual colors?
  facet_wrap("ref") +
  labs(y = 'seconds')
```

![](figures/fit_db_benchmarks-1.png)<!-- -->

``` r
ref_sizes <-
  read_tsv(here('subworkflows/2_fit_reference_db/results/ref_sizes.tsv'))
ref_sizes
```

Potential problem: started with 174k reference sequences in RDP; filter
with the hard parameter leaves only \~6k sequences remaining.

``` r
dataset_sizes <-
  read_tsv(here('subworkflows/1_prep_samples/results/dataset_sizes.tsv'))
dataset_sizes
```

``` r
fractions <- read_tsv(here('subworkflows/2_fit_reference_db/results/fraction_reads_mapped.tsv'))
fractions %>% 
  ggplot(aes(x=dataset, y=fraction_mapped, color=dataset)) +
  geom_boxplot(alpha = 0.5) +
  scale_color_manual(values = dataset_colors) +
  facet_wrap("ref") +
  ylim(0, 1) +
  labs(title="Sequences mapped during closed-reference OptiFit")
```

![](figures/fraction_reads_mapped-1.png)<!-- -->
