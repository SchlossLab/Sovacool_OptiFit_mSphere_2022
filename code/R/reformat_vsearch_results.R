source(snakemake@input[["fcns"]])

log_smk()

sapply(c("bench", "sensspec", "div", "map"), function(x) snakemake@input[[x]]) %>%
  map_dfc(read_tsv) %>%
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
