source(snakemake@input[["fcns"]])

log_smk()

sapply(c("bench", "sensspec", "map", "summary"), function(x) snakemake@input[[x]]) %>%
  map_dfc(read_tsv) %>%
  mutate(
    dataset = snakemake@wildcards[["dataset"]],
    method = snakemake@wildcards[["method"]],
    ref = ifelse(method == "de_novo",
      NA,
      "gg"
    ),
    region = NA,
    tool = "vsearch",
    num_otus = as.double(sobs)
  ) %>%
  write_tsv(snakemake@output[["tsv"]])
