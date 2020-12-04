source(snakemake@input[["fcns"]])
log_smk()
rmarkdown::render(here::here(snakemake@input[["Rmd"]]),
  output_format = snakemake@params[["format"]],
  output_file = here::here(snakemake@output[["file"]])
)
