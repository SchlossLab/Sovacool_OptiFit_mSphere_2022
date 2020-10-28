source(snakemake@input[["fcns"]])
log_smk()

get_seed <- function(x) {
  str_replace(x, ".*seed_(\\d+).*", "\\1")
}

best_list_file <- tibble(
  list_file = snakemake@input[["list"]],
  sensspec_file = snakemake@input[["sensspec"]]
) %>%
  mutate(
    seed = unique(get_seed(list_file), get_seed(sensspec_file)), # should be identical
    sensspec = read_tsv(sensspec_file) %>% pull(sensspec)
  ) %>%
  slice_max(sensspec) %>%
  pull(list_file) %>%
  .[[1]] # in case there's a tie, just pick one

file.copy(best_list_file, snakemake@output[["tsv"]])
