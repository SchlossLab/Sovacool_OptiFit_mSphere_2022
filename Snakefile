" Download references & datasets, process with mothur, and benchmark the OptiFit algorithm "
import yaml

configfile: 'config/config.yaml'
figs_meta_filename = 'config/figures.yaml'
fig_meta = yaml.safe_load(open(figs_meta_filename))

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
        pdf='docs/paper.pdf'

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
        tiff='figures/workflow.tiff'
    params:
        dim=fig_meta['workflow']['dim'].strip('c(').rstrip(')')
    shell:
        """
        dot -T tiff -Gsize={params.dim}\! -Gdpi=300 {input.gv} > {output}
        """

rule plot_results_sum:
    input:
        R='code/R/plot_results_sum.R'
    output:
        tiff='figures/results_sum.tiff'
    params:
        dim=fig_meta['results_sum']['dim']
    script:
        'code/R/plot_results_sum.R'

rule render_pdf:
    input:
        Rmd="paper/paper.Rmd",
        R='code/R/render.R',
        fcns="code/R/functions.R",
        rda=rules.calc_results_stats.output.rda,
        deps=['paper/preamble.tex', 'paper/head.tex',
              'paper/references.bib', 'paper/msystems.csl',
              figs_meta_filename,
              rules.plot_workflow.output,
              rules.plot_results_sum.output
              ]
    output:
        file='docs/paper.pdf'
    params:
        format='pdf_document'
    script:
        'code/R/render.R'
