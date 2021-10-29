Did MCC results change?
================

``` r
library(glue)
library(here)
library(tidyverse)
theme_set(theme_bw())
```

Why does it look like OptiClust performed \~4% worse than OptiFit with
the split strategy?

## Did OptiClust results change?

[March 2021 commit:
`0b80681`](https://github.com/SchlossLab/Sovacool_OptiFit_2021/commits/main/subworkflows/1_prep_samples/results/opticlust_results.tsv)

``` r
opticlust <- c('Nov-2020_5d87338', 'Mar-2021_0b80681') %>% 
  map_dfr(function(file) {
    read_tsv(here('exploratory', '2021', 'did-results-change', 
                glue('opticlust_results_{file}.tsv'))) %>% 
    mutate(file = file)
  })

opticlust %>% group_by(file) %>% summarize(min_mcc = min(mcc),
                                           median_mcc = median(mcc),
                                           max_mcc = max(mcc))
```

    ## # A tibble: 2 × 4
    ##   file             min_mcc median_mcc max_mcc
    ##   <chr>              <dbl>      <dbl>   <dbl>
    ## 1 Mar-2021_0b80681   0.714      0.826   0.846
    ## 2 Nov-2020_5d87338   0.714      0.826   0.845

``` r
opticlust %>% filter(dataset == 'human' | dataset == 'mouse') %>% 
  group_by(file, dataset) %>% 
  summarize(min_mcc = min(mcc),
            median_mcc = median(mcc),
            max_mcc = max(mcc))
```

    ## # A tibble: 4 × 5
    ## # Groups:   file [2]
    ##   file             dataset min_mcc median_mcc max_mcc
    ##   <chr>            <chr>     <dbl>      <dbl>   <dbl>
    ## 1 Mar-2021_0b80681 human     0.818      0.821   0.821
    ## 2 Mar-2021_0b80681 mouse     0.830      0.831   0.831
    ## 3 Nov-2020_5d87338 human     0.818      0.821   0.821
    ## 4 Nov-2020_5d87338 mouse     0.830      0.831   0.831

Nope!

## Did OptiFit results change?

``` r
optifit <- c('Dec-2020_5fcc830',
             'Apr-2021_c6fcdbd', 
             'Jul-2021_4090079') %>% 
  map_dfr(function(file) {
    read_tsv(here('exploratory', '2021', 'did-results-change', 
                glue('optifit_split_results_{file}.tsv'))) %>% 
    mutate(file = file)
  })

optifit %>% filter(dataset == 'human' | dataset == 'mouse') %>% 
  group_by(file, dataset) %>% 
  summarize(min_mcc = min(mcc),
            median_mcc = median(mcc),
            max_mcc = max(mcc))
```

    ## # A tibble: 6 × 5
    ## # Groups:   file [3]
    ##   file             dataset min_mcc median_mcc max_mcc
    ##   <chr>            <chr>     <dbl>      <dbl>   <dbl>
    ## 1 Apr-2021_c6fcdbd human     0.617      0.819   0.867
    ## 2 Apr-2021_c6fcdbd mouse     0.660      0.829   0.856
    ## 3 Dec-2020_5fcc830 human     0.681      0.819   0.860
    ## 4 Dec-2020_5fcc830 mouse     0.697      0.829   0.852
    ## 5 Jul-2021_4090079 human     0.773      0.875   0.899
    ## 6 Jul-2021_4090079 mouse     0.786      0.895   0.911

``` r
optifit %>% 
  filter(dataset == 'human' | dataset == 'mouse') %>% 
  ggplot(aes(x = file, y = mcc, color = dataset)) +
  geom_boxplot()
```

![](figures/optifit-1.png)<!-- -->

Yeah! OptiFit results improved. What changed? I don’t think I’ve updated
mothur since December.
