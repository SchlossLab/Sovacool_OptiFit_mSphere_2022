Nov. 2020

``` r
library(cowplot)
library(ggtext)
library(glue)
library(here)
library(knitr)
library(tidyverse)
theme_set(theme_classic())
color_palette <- RColorBrewer::brewer.pal(4, "Dark2")
dataset_colors <- c(
  human = color_palette[[3]],
  marine = color_palette[[1]],
  mouse = color_palette[[4]],
  soil = color_palette[[2]]
)
set.seed(2018)
```

## *de novo* clustering

``` r
dataset_sizes <- read_tsv(here('subworkflows/1_prep_samples/results/dataset_sizes.tsv'))
opticlust <- read_tsv(here('subworkflows/1_prep_samples/results/opticlust_results.tsv'))
```

## fit to reference databases

``` r
ref_sizes <- read_tsv(here('subworkflows/2_fit_reference_db/results/ref_sizes.tsv'))
optifit_dbs <- read_tsv(here('subworkflows/2_fit_reference_db/results/optifit_dbs_results.tsv')) %>% 
  mutate(mem_mb = max_rss,
         mem_gb = mem_mb / 1024,
         sec = s)
head(optifit_dbs)
```

    ## # A tibble: 6 x 47
    ##   label...1 cutoff numotus     tp      tn     fp     fn sensitivity specificity
    ##       <dbl>  <dbl>   <dbl>  <dbl>   <dbl>  <dbl>  <dbl>       <dbl>       <dbl>
    ## 1      0.03   0.03   35104 4.14e7 3.23e10 6.16e6 1.39e7       0.748        1.00
    ## 2      0.03   0.03   35084 4.22e7 3.23e10 6.37e6 1.31e7       0.763        1.00
    ## 3      0.03   0.03   35116 4.22e7 3.23e10 6.37e6 1.31e7       0.762        1.00
    ## 4      0.03   0.03   35075 4.22e7 3.23e10 6.39e6 1.31e7       0.763        1.00
    ## 5      0.03   0.03   34653 4.21e7 3.23e10 6.33e6 1.33e7       0.76         1.00
    ## 6      0.03   0.03   35119 4.22e7 3.23e10 6.40e6 1.31e7       0.763        1.00
    ## # … with 38 more variables: ppv <dbl>, npv <dbl>, fdr <dbl>, accuracy <dbl>,
    ## #   mcc <dbl>, f1score <dbl>, s <dbl>, `h:m:s` <time>, max_rss <dbl>,
    ## #   max_vms <dbl>, max_uss <dbl>, max_pss <dbl>, io_in <dbl>, io_out <dbl>,
    ## #   mean_load <dbl>, label...25 <dbl>, group <lgl>, nseqs <dbl>, sobs <dbl>,
    ## #   npshannon <dbl>, invsimpson <dbl>, invsimpson_lci <dbl>,
    ## #   invsimpson_hci <dbl>, dataset <chr>, ref <chr>, region <chr>, seed <dbl>,
    ## #   method <chr>, printref <lgl>, fraction_mapped <dbl>, ref_weight <lgl>,
    ## #   ref_frac <lgl>, sample_frac <lgl>, tool <chr>, label...24 <dbl>,
    ## #   mem_mb <dbl>, mem_gb <dbl>, sec <dbl>

``` r
optifit_dbs %>% 
  ggplot(aes(x=method, y=mcc, color=ref)) +
  geom_jitter(alpha = 0.5) + 
  facet_wrap('dataset') +
  ylim(0, 1)
```

![](figures/fit-db_mcc-1.png)<!-- -->

``` r
optifit_dbs %>% filter(method == 'closed') %>% 
  ggplot(aes(x=dataset, y=fraction_mapped, color=ref)) +
  geom_jitter(alpha = 0.5) + 
  ylim(0, 1) +
  labs(title="Sequences mapped during closed-reference OptiFit")
```

![](figures/fit-db_fraction-mapped-1.png)<!-- -->

``` r
optifit_dbs %>% 
  ggplot(aes(x=method, y=sec, color=ref)) +
  geom_jitter(alpha = 0.5) +
  facet_wrap('dataset', scales = 'free')
```

![](figures/fit-db_runtime-1.png)<!-- -->

## fit to subsamples

## vsearch

For reference-based clustering, datasets were fit to the greengenes
database.

``` r
vsearch <- read_tsv(here('subworkflows/4_vsearch/results/vsearch_results.tsv')) %>% 
  mutate(mem_mb = max_rss,
         mem_gb = mem_mb / 1024,
         sec = s)
```

``` r
vsearch %>% 
  pivot_longer(c(mcc, sec), names_to = 'metric') %>% 
  ggplot(aes(x = dataset, y = value, color = method)) +
  geom_boxplot() +
  facet_wrap('metric', scales = 'free', ncol = 1) + 
  #scale_color_manual(values = dataset_colors) +
  theme(axis.title = element_blank())
```

![](figures/vsearch_mcc-1.png)<!-- -->
