
### data prep

set.seed(20210509)
library(cowplot)
library(ggtext)
library(glue)
library(here)
library(knitr)
library(tidyverse)

color_palette <- RColorBrewer::brewer.pal(4, "Dark2")
dataset_colors <- c(
  human = color_palette[[3]],
  marine = color_palette[[1]],
  mouse = color_palette[[4]],
  soil = color_palette[[2]]
)
tri_colors <-  c("darkorange","darkorchid","cyan4") # https://allisonhorst.github.io/palmerpenguins/articles/intro.html
dual_colors <- RColorBrewer::brewer.pal(3, 'Set1')[1:2]
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
select_cols <- function(dat) {
  dat %>% 
    select(dataset, strategy, method, tool, mcc, sec, mem_gb, fraction_mapped)
}
group_sum <- function(dat) {
  dat %>% 
    group_by(dataset, tool, method, strategy) %>% 
    summarize(mcc_median = median(mcc),  # TODO: tidy way to avoid this repetitiveness?
              sec_median = median(sec),
              mem_gb_median = median(mem_gb),
              frac_map_median = median(fraction_mapped))
}

opticlust <- read_tsv(here('subworkflows/1_prep_samples/results/opticlust_results.tsv')) %>% 
  full_join(read_tsv(here('subworkflows/1_prep_samples/results/dataset_sizes.tsv'))) %>% 
  mutate_perf() %>% mutate(fraction_mapped = 1, strategy = method)
sum_opticlust <- opticlust %>% group_sum()
optifit_dbs <- read_tsv(here('subworkflows/2_fit_reference_db/results/optifit_dbs_results.tsv')) %>% 
  mutate_perf()
optifit_split <- read_tsv(here('subworkflows/3_fit_sample_split/results/optifit_split_results.tsv')) %>% 
  mutate_perf()
optifit_all <- list(optifit_dbs %>% 
                      mutate(strategy = glue('database_{ref}')),
                    optifit_split %>% 
                      mutate(strategy = 'self-split')) %>% 
  reduce(full_join)
sum_optifit <- optifit_all %>% group_sum()
vsearch <- read_tsv(here('subworkflows/4_vsearch/results/vsearch_results.tsv')) %>% 
  mutate_perf() %>% 
  mutate(strategy = case_when(
    method == 'de_novo' ~ method,
    TRUE ~ as.character(glue('database_{ref}')))) 
mothur_vsearch <- list(optifit_all, opticlust, vsearch) %>% 
  lapply(select_cols) %>% 
  reduce(bind_rows) %>% 
  mutate(method = as.character(method),
         strategy = as.character(strategy))
sum_all <- mothur_vsearch %>% group_sum()
head(sum_all)

### make the strategy & method labels prettier

sum_all <- sum_all %>%
  mutate(
    strategy = factor(
      case_when(
        strategy == "de_novo"        ~ "_de novo_",
        strategy == 'database_rdp'   ~ "db: RDP",
        strategy == 'database_silva' ~ "db: SILVA",
        strategy == 'database_gg'    ~ "db: Greengenes",
        TRUE                         ~ strategy
      ),
      levels = c('db: RDP', 'db: SILVA', 'db: Greengenes', 
                 'self-split',  '_de novo_')
    ),
    method = factor(
      case_when(method == "de_novo" ~ "_de novo_",
                TRUE                ~ method),
      levels = c('open', 'closed', '_de novo_')
    )
  )

## plot!
plot_quality <- function(dat, y_val, title = '') {
  dat %>% 
    ggplot(aes(strategy, {{ y_val }}, 
               color = tool, 
               shape = method)) + 
    geom_point(size = 3, position = position_dodge(width = 0.4)) +
    facet_wrap(dataset ~ ., nrow=1) +
    scale_shape_manual(values = list(open = 1, closed = 19, `_de novo_` = 17)) +
    scale_color_manual(values = dual_colors) +
    scale_y_continuous(labels = c('0', '0.5', '1'), 
                       breaks = c(0, 0.5, 1),
                       limits = c(0, 1)) +
    coord_flip() +
    labs(x = '', y = '', title = title) +
    theme_bw() +
    theme(legend.position="none",
          axis.text.y = element_markdown())
}

mcc_plot <- sum_all %>% 
  plot_quality(mcc_median, title = "MCC")
frac_plot <- sum_all %>% filter(method == 'closed') %>% 
  plot_quality(frac_map_median, title = "Fraction Mapped")

plot_runtime <- function(dat, yval, title = '') {
  dat %>% 
    ggplot(aes(strategy, {{ yval }}, 
               color = tool, 
               shape = method
    )) + 
    geom_point(size = 3, position = position_dodge(width = 0.4)) +
    facet_wrap(dataset ~ ., nrow = 1, scales = 'free_x') +
    scale_shape_manual(values = list(open = 1, closed = 19, `_de novo_` = 17)) +
    scale_y_log10() +
    scale_color_manual(values = dual_colors) +
    coord_flip() +
    labs(x = '', y = '', title = title) +
    theme_bw() +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "none",
          axis.text.y = element_markdown())
}
runtime_plot <- sum_all %>% plot_runtime(sec_median) + 
  labs(title = 'Runtime (sec)')

shared_legend <- get_legend(mcc_plot +
                              guides(color = guide_legend(nrow = 1)) +
                              theme(legend.position = "bottom",
                                    legend.title = element_blank(),
                                    legend.text = element_markdown())
)

main_plot <- plot_grid(mcc_plot, frac_plot, runtime_plot,
                       ncol = 1, align = 'v', labels = 'AUTO'
) 

plot_grid(main_plot, shared_legend, 
          ncol = 1, rel_heights = c(1, 0.1))

ggsave(here('figures', 'results.tiff'), 
       device = 'tiff', dpi=300,
       width=6, height=6, units='in')

# TODO: reduce panel gaps/margins
# TODO: jitter only needed for some points, e.g. de novo doesn't need jitter. 
# TODO: reverse runtime x scale