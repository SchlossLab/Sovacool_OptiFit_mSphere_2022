source(snakemake@input[["fcns"]])

log_smk()
rmarkdown::render(here::here(snakemake@input[["Rmd"]]))