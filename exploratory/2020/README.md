Exploratory Plots
================
5/27/2020

``` r
library(here)
library(tidyverse)
theme_set(theme_classic())
color_palette <- RColorBrewer::brewer.pal(4, "Dark2")
dataset_colors <- list(human = color_palette[[3]],
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
plot_jitter <- function(.df, ...) {
  .df %>% 
    ggplot(aes(...)) +
    geom_jitter() +
    scale_color_manual(values=dataset_colors)
}
```

``` r
sensspec <- read_tsv(here('subworkflows/2_fit_reference_db/results/sensspec.txt')) %>% 
  relevel_method()
sensspec %>% 
  group_by(dataset, method) %>% 
  plot_jitter(x = method, y = mcc, color = dataset) +
  ylim(0, 1)
```

![](figures/fit_db_sensspec-1.png)<!-- -->

``` r
benchmarks <- read_tsv(here('subworkflows/2_fit_reference_db/results/benchmarks.txt')) %>% 
  relevel_method
benchmarks %>% 
  group_by(dataset, method) %>% 
  plot_jitter(x = method, y = s, color = dataset) +
  labs(y='seconds')
```

![](figures/fit_db_benchmarks-1.png)<!-- -->

``` r
fractions <- read_tsv(here('subworkflows/2_fit_reference_db/results/fraction_reads_mapped.txt'))
fractions %>% 
  plot_jitter(x=dataset, y=fraction_mapped, color=dataset) +
  ylim(0, 1)
```

![](figures/fraction_reads_mapped-1.png)<!-- -->
