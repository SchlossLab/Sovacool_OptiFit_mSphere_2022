source(snakemake@input[["fcns"]])
log_smk()
rmarkdown::render(
  here::here(snakemake@input[["Rmd"]]),
  params = list(snakemake = snakemake),
  output_format = "all",
  output_dir = here::here(snakemake@params[["outdir"]])
)
