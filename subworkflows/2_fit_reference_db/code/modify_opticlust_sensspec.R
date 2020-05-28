sensspec <- readr::read_tsv(snakemake@input[['txt']]) %>%
  dplyr::mutate(dataset = snakemake@params[['dataset']],
         ref = snakemake@params[['ref']],
         region = snakemake@params[['region']],
         seed = snakemake@params[['seed']],
         method = 'de_novo',
         printref = NA,
         iter = NA,
         numotus = NA)
readr::write_tsv(sensspec, snakemake@output[['txt']])