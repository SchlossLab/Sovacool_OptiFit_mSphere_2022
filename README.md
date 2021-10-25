# OptiFit: an improved method for fitting amplicon sequences to existing OTUs

Benchmarking the OptiFit algorithm in [mothur](https://github.com/mothur/mothur).

Find details on how to use OptiFit and descriptions of the parameter
options on the mothur wiki: https://mothur.org/wiki/cluster.fit/.

## The Workflow

The workflow is split into five subworkflows:

- **[0_prep_db](subworkflows/0_prep_db)** — download & preprocess reference
    databases.
- **[1_prep_samples](subworkflows/1_prep_samples)** — download, preprocess, &
    _de novo_ cluster the sample datasets.
- **[2_fit_reference_db](subworkflows/2_fit_reference_db)** — fit datasets to
    reference databases.
- **[3_fit_sample_split](subworkflows/3_fit_sample_split)** — split datasets;
    cluster one fraction _de novo_ and fit the remaining sequences to the
    _de novo_ OTUs.
- **[4_vsearch](subworkflows/4_vsearch)** — run vsearch clustering for
    comparison.

The main workflow ([`Snakefile`](Snakefile)) creates plots from the results of
the subworkflows and renders the [paper](paper).

## Quickstart

1. Before cloning, configure git symlinks:
    ```bash
    git config --global core.symlinks true
    ```
    Otherwise, git will create text files in place of symlinks.
1. Clone this repository.
    ```bash
    git clone https://github.com/SchlossLab/OptiFitAnalysis
    cd OptiFitAnalysis/
    ```
1. Install the dependencies.
    Almost all are listed in the conda environment file.
    ```bash
    conda env create -f config/env.simple.yaml
    conda activate optifit
    ```
    Additionally, I used a custom version of `ggraph` for the algorithm figure.
    You can install it with devtools from R:
    ```r
    devtools::install_github('kelly-sovacool/ggraph', ref = 'iss-297_ggtext')
    ```
    The schtools package is also needed to render the manuscript:
    ```r
    devtools::install_github('SchlossLab/schtools')
    ```
1. Run the entire pipeline.
    Locally:
    ```
    snakemake --cores 4
    ```
    Or on an HPC running slurm:
    ```
    sbatch code/slurm/submit_all.sh
    ```
    (You will first need to edit your email and slurm account info in the
    [submission script](code/slurm/submit.sh)
    and [cluster config](config/cluster.json).)

## Directory Structure

```
.
├── OptiFit.Rproj
├── README.md
├── Snakefile
├── code
│   ├── R
│   ├── bash
│   ├── py
│   ├── slurm
│   └── tests
├── config
│   ├── cluster.json
│   ├── config.yaml
│   ├── config_test.yaml
│   ├── env.export.yaml
│   ├── env.simple.yaml
│   └── slurm
│       └── config.yaml
├── docs
│   ├── paper.md
│   ├── paper.pdf
│   └── slides
├── exploratory
│   ├── 2018_fall_rotation
│   ├── 2019_winter_rotation
│   ├── 2020-05_May-Oct
│   ├── 2020-11_Nov-Dec
│   ├── 2021
│   │   ├── figures
│   │   ├── plots.Rmd
│   │   ├── plots.md
│   ├── AnalysisRoadmap.md
│   └── DeveloperNotes.md
├── figures
├── log
├── paper
│   ├── figures.yaml
│   ├── head.tex
│   ├── msphere.csl
│   ├── paper.Rmd
│   ├── preamble.tex
│   └── references.bib
├── results
│   ├── aggregated.tsv
│   ├── stats.RData
│   └── summarized.tsv
└── subworkflows
    ├── 0_prep_db
    │   ├── README.md
    │   └── Snakefile
    ├── 1_prep_samples
    │   ├── README.md
    │   ├── Snakefile
    │   ├── data
    │   │   ├── human
    │   │       └── SRR_Acc_List.txt
    │   │   ├── marine
    │   │       └── SRR_Acc_List.txt
    │   │   ├── mouse
    │   │       └── SRR_Acc_List.txt
    │   │   └── soil
    │   │       └── SRR_Acc_List.txt
    │   └── results
    │       ├── dataset_sizes.tsv
    │       └── opticlust_results.tsv
    ├── 2_fit_reference_db
    │   ├── README.md
    │   ├── Snakefile
    │   └── results
    │       ├── denovo_dbs.tsv
    │       ├── optifit_dbs_results.tsv
    │       └── ref_sizes.tsv
    ├── 3_fit_sample_split
    │   ├── README.md
    │   ├── Snakefile
    │   └── results
    │       ├── optifit_crit_check.tsv
    │       └── optifit_split_results.tsv
    └── 4_vsearch
        ├── README.md
        ├── Snakefile
        └── results
            ├── vsearch_abbr.md
            └── vsearch_results.tsv

```