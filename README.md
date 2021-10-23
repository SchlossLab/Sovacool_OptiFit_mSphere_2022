# OptiFitAnalysis

Benchmarking the OptiFit algorithm in [mothur](https://github.com/mothur/mothur).

Find more details on how to use OptiFit and descriptions of the parameter
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
    sbatch code/slurm/submit.sh
    ```
    (You will first need to edit your email and slurm account info in the 
    [submission script](code/slurm/submit.sh) 
    and [cluster config](config/cluster.json).)

## Developer Notes

### Project progress

See the [Analysis Roadmap](https://github.com/SchlossLab/OptiFitAnalysis/blob/master/AnalysisRoadmap.md).

### Whitespace in Python

I have my editor set to convert tabs to spaces with a tab length of 4.
In [Atom](https://atom.io) you can use 
[:untabify](https://atom.io/packages/tabs-to-spaces) to convert tabs to spaces.
It's crucial for this to be consistent within Python & Snakemake files.
If you get an error while running snakemake like:
```
Unexpected keyword <word> in rule definition (Snakefile, line <line>)
```
It's likely a whitespace issue.

### Managing software dependencies

I use the [conda](https://conda.io/docs/) package manager to manage software
dependencies.
If you don't already have it, I recommend installing the
[Miniconda](https://conda.io/miniconda.html) Python 3 distribution.
If you're experiencing slowness when solving conda environments, follow the
suggestions [here](https://github.com/bioconda/bioconda-recipes/issues/13774)
and also check out [mamba](https://mamba.readthedocs.io/en/latest/)

#### Create a conda environment

If you plan to run this workflow on the GreatLakes HPC or another 64-bit Linux
machine, you can get an exact replica of my conda environment with:
```
conda env create --file config/env.export.yaml
```

Otherwise, run:
```
conda env create --name optifit --file config/env.simple.yaml
```

#### Activate

Activate the environment before running any code with:
```
conda activate optifit
```

#### Update the environment

Install new packages with:
```
conda install new_package_name
```

Always update the environment file after installing new packages:
```
conda env export > config/env.export.yaml
```

And update the simple version (`config/env.simple.yaml`) with your favorite text
editor. 
The exported file `config/env.export.yaml` contains an exhaustive list of
dependencies and their exact version numbers.
