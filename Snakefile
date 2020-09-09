" Download references & datasets, process with mothur, and benchmark the OptiFit algorithm "

configfile: 'config/config.yaml'

subworkflow prep_db:
    workdir:
        "subworkflows/0_prep_db"

subworkflow prep_samples:
    workdir:
        "subworkflows/1_prep_samples"

subworkflow fit_ref_db:
    workdir:
        "subworkflows/2_fit_reference_db"

rule targets:
    input:
        'paper/paper.pdf'

rule render_paper:
    input:
        Rmd="paper/paper.Rmd",
        pre='paper/header.tex',
        bib='paper/references.bib',
        csl='paper/msystems.csl',
        R='code/R/render.R',
        fcns="code/R/functions.R"
    output:
        pdf='paper/paper.pdf'
    log:
        'log/render_paper.txt'
    script:
        'code/R/render.R'
