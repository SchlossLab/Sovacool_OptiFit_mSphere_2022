rmarkdown::render(
  here::here(snakemake@input[["Rmd"]]),
  params = list(snakemake = snakemake),
  output_format = snakemake@params[["format"]],
  output_file = here::here(snakemake@output[1])
)
