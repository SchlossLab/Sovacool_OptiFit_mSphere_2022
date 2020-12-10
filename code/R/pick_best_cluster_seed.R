source(snakemake@input[["fcns"]])
log_smk()

get_seed <- function(x) {
  str_replace(x, ".*seed_(\\d+).*", "\\1")
}
read_sensspec <- function(filename) {
  read_tsv(filename) %>%
    mutate(seed = get_seed(filename))
}

best_seed <- snakemake@input[["sensspec"]] %>%
  map_dfr(read_sensspec) %>%
  top_n(1, mcc) %>%
  pull(seed) %>%
  .[[1]]

best_list_file <- tibble(list_file = snakemake@input[["list"]]) %>%
  mutate(
    seed = get_seed(list_file)
  ) %>%
  filter(seed == best_seed) %>%
  pull(list_file) %>%
  .[[1]]

file.copy(best_list_file, snakemake@output[["list"]])
