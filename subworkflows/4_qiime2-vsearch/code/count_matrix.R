library(dplyr)

dat <- readr::read_tsv(snakemake@input[["count"]])
seqs <- dat %>% pull(Representative_Sequence)
counts = dat %>% pull(total)
len <- length(seqs)

mat <- matrix(0, nrow = length(seqs), ncol = length(seqs))
rownames(mat) <- seqs
colnames(mat) <- seqs
diag(mat) <- counts

mat %>% 
  as_tibble() %>% 
  mutate(Representative_Sequence = seqs) %>%
  readr::write_tsv(snakemake@output[["count"]])