source(snakemake@input[["fcns"]])
log_smk()
rmarkdown::render(
  here::here(snakemake@input[["Rmd"]]),
  params = list(smk = snakemake),
  output_format = snakemake@params[["format"]],
  output_file = here::here(snakemake@output[["file"]])
)
