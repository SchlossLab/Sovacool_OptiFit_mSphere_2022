rmarkdown::render(
  here::here(snakemake@input[["Rmd"]]),
  params = list(include_figures = snakemake@params[["include_figures"]]),
  output_format = snakemake@params[["format"]],
  output_file = here::here(snakemake@output[1])
)
