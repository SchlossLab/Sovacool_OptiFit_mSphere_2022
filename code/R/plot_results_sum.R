# plot a summary of the results for all datasets & clustering strategies

set.seed(20210509)
library(cowplot)
library(ggtext)
library(glue)
library(here)
library(knitr)
library(tidyverse)
mutate_perf <- function(dat) {
  dat %>%
    mutate(mem_mb = max_rss,
           mem_gb = mem_mb / 1024) %>%
    rename(sec = s)
}
select_cols <- function(dat) {
  dat %>%
    select(dataset, strategy, method, tool, mcc, sec, mem_gb, fraction_mapped)
}

opticlust <- read_tsv(here('subworkflows/1_prep_samples/results/opticlust_results.tsv')) %>%
  full_join(read_tsv(here('subworkflows/1_prep_samples/results/dataset_sizes.tsv'))) %>%
  mutate_perf() %>%
  mutate(strategy = method, fraction_mapped = NA)
optifit_dbs <- read_tsv(here('subworkflows/2_fit_reference_db/results/optifit_dbs_results.tsv')) %>%
  mutate_perf()
optifit_split <- read_tsv(here('subworkflows/3_fit_sample_split/results/optifit_split_results.tsv')) %>%
  filter(ref_frac == 0.5, ref_weight == 'simple') %>%
  mutate_perf()
optifit_all <- list(optifit_dbs %>%
                      mutate(strategy = glue('database_{ref}')),
                    optifit_split %>%
                      mutate(strategy = 'self-split')) %>%
  reduce(full_join)
vsearch <- read_tsv(here('subworkflows/4_vsearch/results/vsearch_results.tsv')) %>%
  mutate_perf() %>%
  mutate(strategy = case_when(
    method == 'de_novo' ~ method,
    TRUE ~ as.character(glue('database_{ref}'))))
mothur_vsearch <- list(optifit_all, opticlust, vsearch) %>%
  lapply(select_cols) %>%
  reduce(bind_rows) %>%
  mutate(method = as.character(method),
         strategy = as.character(strategy)) %>%
  mutate(fraction_mapped = case_when(
    method %>% as.character() != 'closed'  ~ NA_real_,
    TRUE                                   ~ fraction_mapped
  )) %>%
  pivot_longer(c(mcc, fraction_mapped, sec),
               names_to = 'metric') %>%
  mutate(
    metric = factor(
      case_when(
        metric == "mcc"        ~ "MCC",
        metric == 'fraction_mapped'   ~ "Fraction Mapped",
        metric == 'sec'        ~ "Runtime (sec)",
        TRUE                          ~ metric
      ),
      levels = c('MCC', 'Fraction Mapped', 'Runtime (sec)')
    ),
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
    ))

med_iqr <- function(x) {
  return(data.frame(y = median(x),
                    ymin = quantile(x)[2],
                    ymax = quantile(x)[4]))
}

color_list <- list(mothur = RColorBrewer::brewer.pal(3, 'Set1')[1],
                     vsearch = RColorBrewer::brewer.pal(3, 'Set1')[2])
color_labels <- lapply(names(color_list), 
                       function(name) {
                         glue("<span style = 'color:{color_list[[name]]};'>{name}</span>")
                       }
                       ) %>% unlist()

mothur_vsearch  %>%
  ggplot(aes(value, strategy, color = tool, shape = method)) +
  # stat_summary(geom = "linerange",
  #              fun.data = med_iqr,
  #              position = position_dodge(width = 0.4)) +
  stat_summary(geom = 'point',
               fun = median,
               size = 3,
               position = position_dodge(width = 0.4)) +
  facet_grid(dataset ~ metric, scales = 'free', switch = 'x') +
  scale_shape_manual(values = list(open = 1, closed = 19, `_de novo_` = 17)) +
  scale_color_manual(values = color_list,
                     labels = color_labels) +
  labs(x = '', y = '') +
  theme_bw() +
  theme(strip.placement = "outside",
        axis.text.y = element_markdown(),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.position="top",
        legend.margin=margin(t=0, r=0, b=0, l=0, unit='pt'),
        plot.margin=unit(x=c(0,0,0,0),units="pt")
  ) + 
  guides(shape = guide_legend(order = 1),
         colour = guide_legend(override.aes = list(size = -1),
                               order = 2)
  )

dims <- eval(parse(text=snakemake@params[['dim']]))
ggsave(snakemake@output[['tiff']],
       device = 'tiff', dpi=300,
       width=dims[1], height=dims[2], units='in')
