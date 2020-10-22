source(snakemake@input[["fcns"]])

log_smk()

dsrr <- strsplit(snakemake@wildcards[["dataset_ref_region"]], "_") %>% unlist()

sapply(c("bench", "sensspec", "div"), function(x) snakemake@input[[x]]) %>%
  map(read_tsv) %>%
  reduce(full_join) %>%
  mutate(
    dataset = snakemake@wildcards[["dataset"]],
    method = snakemake@wildcards[["method"]],
    ref = ifelse(method == "de_novo",
      NA,
      dsrr[[2]]
    ),
    region = ifelse(method == "de_novo",
      NA,
      paste(dsrr[[3]], dsrr[[4]], sep = "_")
    ),
    tool = "vsearch"
  ) %>%
  write_tsv(snakemake@output[["tsv"]])
