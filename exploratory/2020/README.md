Exploratory Plots
================
5/27/2020

``` r
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
```

# Fitting datasets to a reference database vs *de novo* clustering

``` r
relevel_method <- function(df) {
  df %>%
    mutate(method = fct_relevel(method, c("open", "closed", "de_novo")))
}
```

### Performance as measured by MCC

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

### Performance as measured by runtime

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

### Reference & dataset sizes

``` r
ref_sizes <-
  read_tsv(here('subworkflows/2_fit_reference_db/results/ref_sizes.tsv'))
kable(ref_sizes)
```

| reference | region | num\_seqs | dataset\_filter |
| :-------- | :----- | --------: | :-------------- |
| silva     | v4     |     66531 | human           |
| gg        | v4     |    104943 | human           |
| rdp       | v4     |      6224 | human           |
| silva     | v4     |     66556 | marine          |
| gg        | v4     |    104975 | marine          |
| rdp       | v4     |      6223 | marine          |
| silva     | v4     |     66297 | mouse           |
| gg        | v4     |    104739 | mouse           |
| rdp       | v4     |      6219 | mouse           |
| silva     | v4     |     66593 | soil            |
| gg        | v4     |    104994 | soil            |
| rdp       | v4     |      6224 | soil            |

``` r
dataset_sizes <-
  read_tsv(here('subworkflows/1_prep_samples/results/dataset_sizes.tsv'))
kable(dataset_sizes)
```

| dataset | num\_seqs |
| :------ | --------: |
| human   |    261535 |
| marine  |    161560 |
| mouse   |     68108 |
| soil    |    219752 |

### Fraction of sequences that map to the reference

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
