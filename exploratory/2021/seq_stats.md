2021-10-22

# Similarity & Abundance of Sequences

Each sequence has an absolute abundance and similarity value, where
similarity is the count of other sequences within the dataset that it is
similar to per the 0.03 distance threshold. I was wondering if the
distributions of those values would have some pattern that corresponded
to what we see with the MCC scores for weighting by those methods. There
isn’t anything obviously weird or wrong. For all datasets most sequences
are very low abundance (expected). Both marine and soil have a much
lower median similarity value than human and mouse, and the marine and
soil datasets also have a huge spread in MCC between open & closed with
similarity weighting.

``` r
set.seed(20211021)
library(glue)
library(here)
library(tidyverse)
theme_set(theme_bw())
```

``` r
seq_stats <- c("human","marine","mouse","soil") %>% map_dfr(function(x) {
  read_tsv(glue('subworkflows/1_prep_samples/data/{x}/seq_stats.tsv')) %>% 
    mutate(dataset = x)
}) %>% pivot_longer(c(similarity, abundance), names_to = 'statistic')
```

### histograms

``` r
seq_stats %>% 
  ggplot(aes(x=value, fill=statistic)) +
  geom_histogram() +
  facet_grid(dataset ~ statistic, scales='free')
```

![](figures/hist_sim_abun-1.png)<!-- -->

``` r
seq_stats %>% 
  ggplot(aes(x=value, fill=statistic)) +
  geom_histogram() +
  scale_x_log10() +
  facet_wrap("dataset")
```

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning: Removed 34675 rows containing non-finite values (stat_bin).

![](figures/hist_sim_abun_log10-1.png)<!-- -->

### boxplots

``` r
seq_stats %>% 
  ggplot(aes(x=value, y=dataset, fill=statistic)) +
  geom_boxplot() +
  facet_wrap('statistic', scales = 'free')
```

![](figures/box_sim_abun-1.png)<!-- -->

``` r
seq_stats %>% 
  ggplot(aes(x=value, y=dataset, fill=statistic)) +
  geom_boxplot() +
  scale_x_log10() +
  facet_wrap('statistic', scales = 'free')
```

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning: Removed 34675 rows containing non-finite values (stat_boxplot).

![](figures/box_sim_abun_log10-1.png)<!-- -->

### abundance vs similarity

``` r
seq_stats %>% 
  pivot_wider(names_from = statistic) %>% 
  ggplot(aes(x = abundance, y = similarity, color = dataset)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  facet_wrap('dataset')
```

![](figures/point_abun_sim-1.png)<!-- -->

## relative abundance & simlarity

``` r
totals <- seq_stats %>%
  pivot_wider(names_from = statistic) %>% 
  group_by(dataset) %>% 
  summarize(tot_simi = sum(similarity),
            tot_abun = sum(abundance))

rel_stats <- c("human", "marine", "mouse", "soil") %>% map_dfr(function(x) {
  read_tsv(glue('subworkflows/1_prep_samples/data/{x}/seq_stats.tsv')) %>% 
    mutate(dataset = x,
           rel_sim = similarity / (totals %>% 
                                     filter(dataset == x) %>% 
                                     pull(tot_simi)),
           rel_abun = abundance / (totals %>% 
                                     filter(dataset == x) %>% 
                                     pull(tot_abun))
           )
})
```

### rel abun vs rel sim

``` r
rel_stats %>% 
  ggplot(aes(x = rel_abun, y = rel_sim, color = dataset)) +
  geom_point(alpha=0.7) +
  facet_wrap('dataset')
```

![](figures/point_rel_abun_sim-1.png)<!-- -->

``` r
rel_stats %>% 
  ggplot(aes(x = rel_abun, y = rel_sim, color = dataset)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  facet_wrap('dataset')
```

![](figures/point_rel_abun_sim_log10-1.png)<!-- -->

### relative boxplots

``` r
rel_stats %>% pivot_longer(c(rel_sim, rel_abun), names_to = 'statistic') %>% 
  ggplot(aes(x=value, y=dataset, fill=statistic)) +
  geom_boxplot() +
  facet_wrap('statistic', scales = 'free')
```

![](figures/box_rel_sim_abun-1.png)<!-- -->

``` r
rel_stats %>% pivot_longer(c(rel_sim, rel_abun), names_to = 'statistic') %>% 
  ggplot(aes(x=value, y=dataset, fill=statistic)) +
  geom_boxplot() +
  scale_x_log10() +
  facet_wrap('statistic', scales = 'free')
```

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning: Removed 34675 rows containing non-finite values (stat_boxplot).

![](figures/box_rel_sim_abun_log10-1.png)<!-- -->

### relative histograms

``` r
rel_stats %>%  pivot_longer(c(rel_sim, rel_abun), names_to = 'statistic') %>% 
  ggplot(aes(x=value, fill=statistic)) +
  geom_histogram() +
  facet_wrap("dataset")
```

![](figures/hist_rel_sim_abun-1.png)<!-- -->

``` r
rel_stats %>%  pivot_longer(c(rel_sim, rel_abun), names_to = 'statistic') %>% 
  ggplot(aes(x=value, fill=statistic)) +
  geom_histogram() +
  scale_x_log10() +
  facet_wrap("dataset")
```

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning: Removed 34675 rows containing non-finite values (stat_bin).

![](figures/hist_rel_sim_abun_log10-1.png)<!-- -->

## summary stats

``` r
rel_stats %>% 
  group_by(dataset) %>% 
  summarize_at(vars(abundance), 
               list(min=min, median=median, mean=mean, max=max)) %>% 
  knitr::kable()
```

| dataset | min | median |      mean |     max |
|:--------|----:|-------:|----------:|--------:|
| human   |   1 |      1 | 158.58380 | 1286715 |
| marine  |   1 |      1 |  15.16795 |  156918 |
| mouse   |   1 |      1 |  78.77425 |  264594 |
| soil    |   1 |      1 |   5.71222 |   42712 |

``` r
rel_stats %>% 
  group_by(dataset) %>% 
  summarize_at(vars(similarity), 
               list(min=min, median=median, mean=mean, max=max)) %>% 
  knitr::kable()
```

| dataset | min | median |     mean |  max |
|:--------|----:|-------:|---------:|-----:|
| human   |   0 |    239 | 614.9974 | 4159 |
| marine  |   0 |     26 | 374.8955 | 5495 |
| mouse   |   0 |    125 | 380.7253 | 1941 |
| soil    |   0 |     22 | 136.0720 | 3778 |
