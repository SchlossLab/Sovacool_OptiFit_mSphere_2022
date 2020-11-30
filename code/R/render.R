source(snakemake@input[["fcns"]])

log_smk()
rmarkdown::render(here::here(snakemake@input[["Rmd"]]),
                  output_format = 'pdf_document',
                  output_file = here::here('docs', 'paper.pdf'))
rmarkdown::render(here::here(snakemake@input[["Rmd"]]), 
                  output_format = 'distill::distill_article', 
                  output_file = here::here('docs','index.html'))
