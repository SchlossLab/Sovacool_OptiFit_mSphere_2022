source(snakemake@input[["fcns"]])

log_smk()

sapply(c("bench", "sensspec", "div", "map"), function(x) snakemake@input[[x]]) %>%
  map(read_tsv) %>%
  reduce(bind_cols) %>%
  mutate(
    dataset = snakemake@wildcards[["dataset"]],
    method = snakemake@wildcards[["method"]],
    ref = ifelse(method == "de_novo",
      NA,
      "gg"
    ),
    region = NA,
    tool = "vsearch"
  ) %>%
  write_tsv(snakemake@output[["tsv"]])
