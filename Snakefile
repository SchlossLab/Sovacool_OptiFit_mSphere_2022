" Download references & datasets, process with mothur, and benchmark the OptiFit algorithm "

configfile: 'config/config.yaml'

subworkflow prep_db:
    workdir:
        "subworkflows/0_prep_db/"

subworkflow prep_samples:
    workdir:
        "subworkflows/1_prep_samples/"

subworkflow fit_ref_db:
    workdir:
        "subworkflows/2_fit_reference_db/"

subworkflow fit_subset:
    workdir:
        "subworkflows/3_fit_sample_subset/"

subworkflow vsearch:
    workdir:
        "subworkflows/4_vsearch/"

rule targets:
    input:
        'paper/paper.pdf',
        prep_samples('results/opticlust_results.tsv'),
        fit_ref_db('results/optifit_dbs_results.tsv'),
        fit_subset('results/optifit_subset_results.tsv'),
        vsearch('results/vsearch_results.tsv')

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
