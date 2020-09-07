source(snakemake@input[['fcns']])
merge_results("sensspec")
merge_results("bench")
merge_results('mapped')