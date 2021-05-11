Comparing OTU quality across clustering strategies
================
2021-05-11

### data prep

``` r
set.seed(20210509)
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
tri_colors <-  c("darkorange","darkorchid","cyan4") # https://allisonhorst.github.io/palmerpenguins/articles/intro.html
mutate_perf <- function(dat) {
  dat %>% 
    mutate(mem_mb = max_rss,
           mem_gb = mem_mb / 1024) %>% 
    rename(sec = s)
}
plot_denovo_hline <- function(yint, dat = sum_opticlust) {
  list(geom_hline(data = dat, aes(yintercept = {{yint}})),
       labs(caption = "Black line: _de novo_ clustering"),
       theme(plot.caption = element_markdown())
  )
}
```

``` r
opticlust <- read_tsv(here('subworkflows/1_prep_samples/results/opticlust_results.tsv')) %>% 
  full_join(read_tsv(here('subworkflows/1_prep_samples/results/dataset_sizes.tsv'))) %>% 
  mutate_perf()
sum_opticlust <- opticlust %>% 
  group_by(dataset) %>% 
  summarize(mcc_median = median(mcc),
            sec_median = median(sec),
            mem_gb_median = median(mem_gb)) %>% 
  mutate(frac_map_median = 1)

optifit_dbs <- read_tsv(here('subworkflows/2_fit_reference_db/results/optifit_dbs_results.tsv')) %>% 
  mutate_perf()
optifit_split <- read_tsv(here('subworkflows/3_fit_sample_split/results/optifit_split_results.tsv')) %>% 
  mutate_perf()
optifit_all <- list(optifit_dbs %>% 
                   mutate(strategy = glue('database_{ref}')),
                 optifit_split %>% 
                   mutate(strategy = 'self-split')) %>% 
  reduce(full_join)

sum_optifit <- optifit_all %>% 
  group_by(dataset, strategy, method) %>% 
  summarize(n = n(),
            mcc_median = median(mcc),  # TODO: tidy way to avoid this repetitiveness?
            sec_median = median(sec),
            mem_gb_median = median(mem_gb),
            frac_map_median = median(fraction_mapped))
head(sum_optifit)
```

    ## # A tibble: 6 x 8
    ## # Groups:   dataset, strategy [3]
    ##   dataset strategy       method     n mcc_median sec_median mem_gb_median
    ##   <chr>   <glue>         <chr>  <int>      <dbl>      <dbl>         <dbl>
    ## 1 human   database_gg    closed   100      0.800       606.          5.38
    ## 2 human   database_gg    open     100      0.815       899.         20.3 
    ## 3 human   database_rdp   closed   100      0.597       476.          5.10
    ## 4 human   database_rdp   open     100      0.819       991.         20.0 
    ## 5 human   database_silva closed   100      0.780       549.          5.22
    ## 6 human   database_silva open     100      0.817       886.         20.1 
    ## # … with 1 more variable: frac_map_median <dbl>

# using facet\_grid()

``` r
sum_optifit %>% 
  pivot_longer(c(frac_map_median, mcc_median), names_to = 'quality_metric') %>% 
  ggplot(aes(strategy, value, color = method)) + 
  plot_denovo_hline(yint = value, dat = sum_opticlust %>% 
                      pivot_longer(c(frac_map_median, mcc_median), 
                                   names_to = 'quality_metric')) +
  geom_point() +
  facet_grid(quality_metric ~ dataset) +
  scale_color_manual(values = tri_colors) +
  coord_flip() +
  labs(x = '', y = '')  + 
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom")
```

![](figures/otu-quality_facet-1.png)<!-- -->

Problems:

-   The sideways grid labels are hard to read.
-   There’s too much whitespace around the legend.

## using cowplot

``` r
plot_quality <- function(dat, y_val, title = '') {
  dat %>% 
    ggplot(aes(strategy, {{ y_val }}, color = method)) + 
    geom_hline(data = sum_opticlust, aes(yintercept = {{ y_val }})) +
    geom_point() +
    facet_wrap(dataset ~ ., nrow=1) +
    scale_color_manual(values = tri_colors) +
    ylim(0, 1) +
    coord_flip() +
    labs(x = '', y = '', title = title) + 
    theme(legend.position="none")
}
mcc_plot <- sum_optifit %>% 
  plot_quality(mcc_median, title = 'Median MCC')
frac_plot <- sum_optifit %>% 
  plot_quality(frac_map_median, title = 'Median fraction mapped')  + 
  labs(caption = "Black line: _de novo_ clustering") +
  theme(plot.caption = element_markdown(hjust = 0),
        plot.caption.position = 'plot')

shared_legend <- get_legend(mcc_plot + 
                              guides(color = guide_legend(nrow = 1)) +
                              theme(legend.position = "bottom")
                            )

main_plot <- plot_grid(mcc_plot, frac_plot, 
                       ncol = 1, align = 'v', labels = 'AUTO')

plot_grid(main_plot, shared_legend, 
          ncol = 1, rel_heights = c(1, 0.1))
```

![](figures/otu-quality_cowplot-1.png)<!-- -->

TODO:

-   Get caption & legend side-by-side. Maybe use `ggannotate` instead of
    `labs(caption)` to accomplish this?
-   Remove excess whitespace around the legend.
    <https://wilkelab.org/cowplot/articles/shared_legends.html>
-   Adjust figure dimensions.

### Feedback from lab meeting

-   [x] vertical jitter
-   [x] frac map only makes sense for closed
-   [x] dotted line for *de novo*
-   [x] filled/unfilled for closed/open
-   [x] box around facets
-   [x] bigger points
-   [ ] set theme elements to add space between facets:
    `+ theme(panel.spacing = unit(c(3,3),'lines'))`
-   generally: build for the dimensions required by the journal

``` r
plot_quality <- function(dat, y_val, title = '') {
  dat %>% 
    ggplot(aes(strategy, {{ y_val }}, 
               #color = method, 
               shape = method)) + 
    facet_wrap(dataset ~ ., nrow=1) +
    scale_shape_manual(values = list(open = 1, closed = 19)) +
    #scale_color_manual(values = tri_colors) +
    scale_y_continuous(labels = c('0', '0.5', '1'), 
                       breaks = c(0, 0.5, 1),
                       limits = c(0, 1)) +
    coord_flip() +
    labs(x = '', y = '', title = title) +
    theme_bw()
}
mcc_plot <- sum_optifit %>% 
  plot_quality(mcc_median, title = 'Median MCC') +
    geom_point(size = 3,
               position = position_dodge(width = 0.2)) +
    geom_hline(data = sum_opticlust, 
               aes(yintercept = mcc_median), 
               linetype = 'dashed', color = 'red')
frac_plot <- sum_optifit %>% filter(method == 'closed') %>% 
  plot_quality(frac_map_median, title = 'Median fraction mapped')  + 
    geom_point(size = 3) +
  labs(caption = "dashed line: _de novo_ clustering") +
  theme(plot.caption = element_markdown(hjust = 0),
        plot.caption.position = 'plot')

shared_legend <- get_legend(mcc_plot + 
                              guides(color = guide_legend(nrow = 1)) +
                              theme(legend.position = "bottom")
                            )

main_plot <- plot_grid(mcc_plot + 
    theme(legend.position="none"), frac_plot + 
    theme(legend.position="none"), 
                       ncol = 1, align = 'v', labels = 'AUTO')

plot_grid(main_plot, shared_legend, 
          ncol = 1, rel_heights = c(1, 0.1))
```

![](figures/otu-quality_cowplot-revised-1.png)<!-- -->
