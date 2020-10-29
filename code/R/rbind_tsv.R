source(snakemake@input[["fcns"]])

log_smk()
merge_results("tsv")
