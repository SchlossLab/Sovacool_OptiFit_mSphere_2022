" Download references & datasets, process with mothur, and benchmark the OptiFit algorithm "
import os
import yaml

configfile: 'config/config.yaml'
figs_meta_filename = 'paper/figures.yaml'
fig_meta = yaml.safe_load(open(figs_meta_filename))
img_dpi = 300
workflow_dim = fig_meta['workflow']['dim'].strip('c(').rstrip(')')
workflow_width = float(workflow_dim.split(',')[0]) * img_dpi
workflow_height = float(workflow_dim.split(',')[1]) * img_dpi

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

rule paper:
    input:
        pdf='docs/paper.pdf',
        md='paper/paper.md'

rule subtargets: # it takes a long time to build the DAG for some of these
    input:
        opticlust=prep_samples('results/opticlust_results.tsv'),
        optifit_db=fit_ref_db('results/optifit_dbs_results.tsv'),
        optifit_split=fit_split('results/optifit_split_results.tsv'),
        vsearch=vsearch('results/vsearch_results.tsv')

rule summarize_results:
    input:
        R='code/R/summarize_results.R',
        opticlust='subworkflows/1_prep_samples/results/opticlust_results.tsv',
        optifit_db='subworkflows/2_fit_reference_db/results/optifit_dbs_results.tsv',
        optifit_split='subworkflows/3_fit_sample_split/results/optifit_split_results.tsv',
        vsearch='subworkflows/4_vsearch/results/vsearch_results.tsv'
    output:
        agg='results/aggregated.tsv',
        sum='results/summarized.tsv',
        vsearch='subworkflows/4_vsearch/results/vsearch_abbr.md'
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

rule plot_algorithm:
    input:
        R='code/R/plot_algorithm_diagram.R'
    output:
        tiff='figures/algorithm.tiff'
    params:
        dim=fig_meta['algorithm']['dim']
    script:
        'code/R/plot_algorithm_diagram.R'

rule plot_workflow: # https://stackoverflow.com/a/20536144/5787827
    input:
        gv='figures/workflow.gv'
    output:
        tmp=temp('figures/workflow.tmp.tiff'),
        tiff='figures/workflow.tiff'
    params:
        dim=workflow_dim,
        width=workflow_width,
        height=workflow_height
    shell:
        """
        dot -T tiff -Gsize={params.dim}\! -Gdpi=300 {input.gv} > {output.tmp}
        convert {output.tmp} -gravity center -background white -extent {params.width}x{params.height} {output.tiff}
        """

rule plot_results_sum:
    input:
        R='code/R/plot_results_sum.R',
        dat=rules.summarize_results.output
    output:
        tiff='figures/results_sum.tiff'
    params:
        dim=fig_meta['results_sum']['dim']
    script:
        'code/R/plot_results_sum.R'

rule plot_results_split:
    input:
        R='code/R/plot_results_split.R',
        dat=rules.summarize_results.output
    output:
        tiff='figures/results_split.tiff'
    params:
        dim=fig_meta['results_split']['dim']
    script:
        'code/R/plot_results_split.R'

rule render_markdown:
    input:
        Rmd="paper/paper.Rmd",
        R='code/R/render.R',
        rda=rules.calc_results_stats.output.rda,
        deps=['paper/preamble.tex', 'paper/head.tex',
              'paper/references.bib', 'paper/msphere.csl',
              figs_meta_filename,
              rules.plot_algorithm.output,
              rules.plot_workflow.output,
              rules.plot_results_sum.output,
              rules.plot_results_split.output
              ]
    output:
        md='paper/paper.md'
    params:
        format='github_document'
    script:
        'code/R/render.R'

rule render_pdf:
    input:
        Rmd="paper/paper.Rmd",
        R='code/R/render.R',
        md_output=rules.render_markdown.output.md
    output:
        pdf='docs/paper.pdf'
    params:
        format='pdf_document'
    script:
        'code/R/render.R'

rule create_test_data:
    input:
        py='code/py/create_test_uc_files.py'
    output:
        'code/tests/data/closed.uc',
        'code/tests/data/denovo.uc',
        'code/tests/data/oracle_open.list'
    script:
        'code/py/create_test_uc_files.py'

rule test_R_code:
    input:
        R='code/tests/testthat.R',
        scripts=[os.path.join('code/R', file) for file in os.listdir('code/R')]
    script:
        'code/tests/testthat.R'

rule test_Python_code:
    input:
        py='code/tests/test_python.py',
        scripts=[os.path.join('code/py', file) for file in os.listdir('code/py')],
        dat=rules.create_test_data.output
    shell:
        'python -m code.tests.test_python'

onsuccess:
    print("🎉 workflow complete!")

onerror:
    print("⛔️ something went wrong...")
