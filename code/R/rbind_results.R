source(snakemake@input[["fcns"]])

log_smk()
merge_results("benchmarks")
merge_results("sensspec")
merge_results("div")