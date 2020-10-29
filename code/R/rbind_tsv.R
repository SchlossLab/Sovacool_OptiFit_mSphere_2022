source(snakemake@input[["fcns"]])

log_smk()
rbind_all("tsv")
