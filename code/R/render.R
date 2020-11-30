source(snakemake@input[["fcns"]])

log_smk()
rmarkdown::render(here::here(snakemake@input[["Rmd"]]),
                  output_format = 'pdf_document')
rmarkdown::render(here::here(snakemake@input[["Rmd"]]), 
                  output_format = 'distill::distill_article', 
                  output_dir = 'docs/')