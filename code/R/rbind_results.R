source(snakemake@input[["fcns"]])

log_smk()
rbind_all("benchmarks")
rbind_all("sensspec")
rbind_all("div")
