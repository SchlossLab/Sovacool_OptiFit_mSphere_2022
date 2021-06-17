" Download references & datasets, process with mothur, and benchmark the OptiFit algorithm "

configfile: 'config/config.yaml'

subworkflow prep_db:
    workdir:
        "subworkflows/0_prep_db/"
    configfile:
        config['configpath']

subworkflow prep_samples:
    workdir:
        "subworkflows/1_prep_samples/"
    configfile:
        config['configpath']

subworkflow fit_ref_db:
    workdir:
        "subworkflows/2_fit_reference_db/"
    configfile:
        config['configpath']

subworkflow fit_split:
    workdir:
        "subworkflows/3_fit_sample_split/"
    configfile:
        config['configpath']

subworkflow vsearch:
    workdir:
        "subworkflows/4_vsearch/"
    configfile:
        config['configpath']

rule render_paper:
    input:
        pdf='docs/paper.pdf',
        html='docs/paper.html'

rule subtargets:
    input:
        opticlust=prep_samples('results/opticlust_results.tsv'),
        optifit_db=fit_ref_db('results/optifit_dbs_results.tsv'),
        optifit_split=fit_split('results/optifit_split_results.tsv'),
        vsearch=vsearch('results/vsearch_results.tsv')

rule summarize_results:
    input:
        R='code/R/summarize_results.R'
        #opticlust=prep_samples('results/opticlust_results.tsv'),
        #optifit_db=fit_ref_db('results/optifit_dbs_results.tsv'),
        #optifit_split=fit_split('results/optifit_split_results.tsv'),
        #vsearch=vsearch('results/vsearch_results.tsv')
    output:
        agg='results/aggregated.tsv',
        sum='results/summarized.tsv'
    script:
        'code/R/summarize_results.R'

rule calc_results_stats:
    input:
        R='code/R/calc_results_stats.R',
        tsv_sum=rules.summarize_results.output.sum,
        tsv2='subworkflows/2_fit_reference_db/results/denovo_dbs.tsv'
    output:
        rda="results/stats.RData"
    script:
        'code/R/calc_results_stats.R'

rule plot_workflow:
    input:
        gv='figures/workflow.gv'
    output:
        'figures/workflow.tiff'
    shell:
        """
        dot -T tiff -Gsize=4,3\! -Gdpi=300 {input.gv} > {output}
        """

deps = ['paper/references.bib',
        'paper/msystems.csl',
        rules.calc_results_stats.output.rda,
        rules.plot_workflow.output]

rule render_pdf:
    input:
        Rmd="paper/paper.Rmd",
        pre=['paper/preamble.tex', 'paper/head.tex', 'paper/tail.tex'],
        R='code/R/render.R',
        fcns="code/R/functions.R",
        deps=deps
    output:
        file='docs/paper.pdf'
    params:
        format='pdf_document'
    script:
        'code/R/render.R'

rule render_html:
    input:
        Rmd="paper/paper.Rmd",
        R='code/R/render.R',
        fcns="code/R/functions.R",
        deps=deps
    output:
        file='docs/paper.html'
    params:
        format="distill::distill_article"
    script:
        'code/R/render.R'

