# OptiFitAnalysis

## Quickstart

1. Install `mothur` version `1.42.0` and all the dependencies listed in `config/environment.yaml`.
2. Clone this repository.
    ```
    git clone https://github.com/SchlossLab/OptiFitAnalysis
    cd OptiFitAnalysis/
    ```
3. Edit the [configuration](https://github.com/SchlossLab/OptiFitAnalysis#snakemake-configuration) file as needed with your preferred text editor.
    ```
    vi config/config.yaml
    ```
4. [Run](https://github.com/SchlossLab/OptiFitAnalysis#running-a-workflow) the pipeline.

    Locally:
    ```
    snakemake
    ```
    Or on the cluster:
    ```
    qsub code/pbs_scripts/optifit-analysis.pbs
    ```

5. Take a look at the results.
    * Each dataset has its own directory in `results/dataset-as-reference` and `results/silva-as-reference`. Within those directories:
        * The main output file is aggregated in `{dataset}/aggregate.sensspec`.
        * Results figures are in `{dataset}/figures/`.

6. Visualize the workflow.
    ```
    snakemake --dag | dot -Tsvg > results/workflows/dag.svg
    ```
    You can open the svg file in your GUI web browser to view it.


## Managing software dependencies

I'm using the [conda](https://conda.io/docs/) package manager to manage dependencies for this project.
If you don't already have it, I recommend installing the [Miniconda](https://conda.io/miniconda.html) Python 3 distribution.
[Here's a link](https://conda.io/docs/_downloads/conda-cheatsheet.pdf) to a helpful cheatsheet for using conda.

### Create a conda environment

If you plan to run this workflow on Flux or other 64-bit Linux machine,
you can get an exact replica of my environment with:
```
conda create --name OptiFitAnalysis --file config/environment.export.txt
```

Otherwise, run:
```
conda env create --name OptiFitAnalysis --file config/environment.yaml
```

The conda project is [in the process of merging conda-env into conda](https://groups.google.com/a/continuum.io/forum/#!topic/conda/EBVVtS8bNRA),
so the syntax for these commands may change.

## Activate

Activate the environment before running any code with:
```
source activate OptiFitAnalysis
```
Be sure to activate the environment from the login node before submitting jobs on the cluster.

## Update

Install new packages with:
```
conda install new_package_name
```

Periodically update the dependencies list:
```
conda list --export > config/environment.export.txt
```
(We update the less-detailed `config/environment.yaml` manually when we install a new package.)

Almost all dependencies are listed in `config/environment.export.txt`, which will be installed by conda when you create the environment.
The exception to this is the `mothur` program -- we're using a different versoin than what's currently available in the bioconda channel.
If you're a member of the Schloss Lab and you're running this analysis on Flux, you can use the mothur binary here: `/nfs/turbo/schloss-lab/bin/mothur-1.42.0/mothur`.
Otherwise, [download the precompiled binary](https://github.com/mothur/mothur/releases) or
[compile from source](https://github.com/mothur/mothur/blob/master/INSTALL.md).
Be sure to use `mothur` version `1.42.0` or higher.

## Snakemake Configuration

The Snakemake workflow relies on a configuration file in `YAML` format. Our default is `config/config.yaml`:

- Set `mothur_bin` to the path to your mothur binary if you're not using the default one on Flux.
- To run the workflow with just a subset of the input data for debugging purposes:
    - Set `subsample_test` to `True`.
    - Set `subsample_size` to the number of sequences you want to use.
- Set `dataset-as-reference` to `True` to run the workflow in `code/analysis/optifit-dataset-as-ref.smk`.
- Set `silva-as-reference` to `True` to run the workflow in `code/analysis/optifit-silva-ref.smk`.
    - Update the version number in `silva_db_version` if needed.
- Add or remove samples in the `datasets` list as needed.
- If using `dataset-as-reference`:
    - `weights` are the methods used to subsample the dataset.
    - `methods` are for controlling `closed`- or `open`-reference clustering in optifit.
    - `printrefs` must be `t` or `f` (mothur parameter) for including the reference data in the output files from optifit.
    - `reference_fractions` controls the different sizes to subsample the dataset to.
    - `iterations` is the number of times the dataset will be randomly subsampled for each fraction.
    - `replicates` is the number of times OptiClust and OptiFit will run on each sample.


## Snakemake Workflows

### Running a workflow

If your workflow is a file named `Snakefile` in your current working dir, run it with:
```
snakemake
```

Otherwise, specify the path to the snakemake file:
```
snakemake -s /path/to/snakefile
```

Override any default `configfile` specified in the workflow with:
```
snakemake --configfile /path/to/config.yaml
```

#### Cores

Run a snakemake workflow with 2 cores:
```
snakemake -j 2
```
Any independent rules are then run in parallel. Without the `-j` or `--cores` flag, snakemake defaults to using only 1 core.
Run a workflow with as many cores as are available with:
```
snakemake -j
```

#### Force run

Snakemake will only run jobs for rules whose output files do not exist and haven't been modified since the last run.
To override this behavior, force a specific rule to run:
```
snakemake --forcerun rule_name
```

Or force all rules to run:
```
snakemake --forceall
```

#### Dry run

Do a dry run to see which jobs snakemake would run without actually running them:
```
snakemake --dryrun
```
Before committing changes or submitting jobs to the cluster, test your snakefile for syntax errors with a dry run.

#### On the cluster

Edit the cluster config file `cluster.json` with your email for PBS to send job notifications.
On the cluster, create a PBS script and append the following command:
```
snakemake --profile pbs-torque
```
Then submit the PBS script to the cluster with `qsub`. See `code/pbs_scripts` for examples.

Profiles for other cluster systems are available in the [snakemake profiles GitHub](https://github.com/snakemake-profiles/doc) such as [Slurm](https://github.com/Snakemake-Profiles/slurm). You may need to edit `cluster.json` for use with a different cluster system profile.

### Visualizing the DAG

Snakemake creates an image representing the directed acyclic graph (DAG) for a workflow with the following command:
```
snakemake --dag | dot -Tsvg > results/workflows/dag.svg
```

Here's a small example DAG:

![optifit-pipeline.config_test.dag.svg](https://github.com/SchlossLab/OptiFitAnalysis/blob/master/results/workflows/optifit-pipeline.config_test.dag.svg)

## Developer Notes

### Whitespace in Python

I have my editor set to convert tabs to spaces with a tab length of 4.
In [Atom](https://atom.io) you can use [:untabify](https://atom.io/packages/tabs-to-spaces) to convert tabs to spaces.
It's crucial for this to be consistent within Python & Snakemake files.
If you get an error while running snakemake like:
```
Unexpected keyword <word> in rule definition (Snakefile, line <line>)
```
It's likely a whitespace issue.

### Project progress

See the [Analysis Roadmap](https://github.com/SchlossLab/OptiFitAnalysis/blob/master/AnalysisRoadmap.md).

### Directory structure

```
.
├── AnalysisRoadmap.md                   # Project progress
├── OptiFitAnalysis.Rproj
├── README.md
├── Snakefile                            # The master snakemake workflow -- ties everything together
├── benchmarks                           # Runtime & mem benchmarks generated by snakemake
│   └── full_dag_usage.txt
├── cluster.json
├── code
│   ├── analysis
│   │   ├── optifit-dataset-as-ref.smk   # Workflow to use the dataset as its own reference
│   │   ├── plot_sensspec.R
│   │   └── weighted_subsample.R
│   ├── data_processing                  # Code for downloading & processing raw data
│   │   ├── get-references.smk
│   │   ├── human.R
│   │   ├── human.batch
│   │   ├── marine.R
│   │   ├── marine.batch
│   │   ├── mice.R
│   │   ├── mice.batch
│   │   ├── soil.R
│   │   ├── soil.batch
│   │   └── testset-subsample.smk
│   ├── legacy                           # Code no longer in use from before Snakemake
│   │   ├── opti_pipeline.sh
│   │   ├── opti_pipeline_flux.sh
│   │   ├── optifit_test.sh
│   │   ├── plot_sensspec_flux.sh
│   │   ├── references.sh
│   │   └── testset-subsample.sh
│   └── pbs_scripts                      # For running jobs on the Flux cluster
│       ├── config_test.pbs
│       ├── optifit-analysis.pbs         # Default script to submit
│       └── tests
│           ├── 1job.pbs
│           ├── mothur-test.pbs
├── config
│   ├── cluster_mumphrey.json
│   ├── cluster_sovacool.json
│   ├── config.yaml                      # Default config file for Snakemake
│   ├── config_test.yaml
│   ├── environment.export.txt           # Exhaustive list of dependencies for use on Flux
│   └── environment.yaml                 # Minimum list of dependencies for other computing environments
├── data                                 # Raw and processed data
│   ├── human
│   ├── marine
│   ├── mice
│   ├── references
│   └── soil
├── exploratory                          # Exploratory analyses & results
│   └── fall_rotation_2018
│       ├── MouseSensspecFigures_2018-12-19.Rmd
│       ├── debug
│       │   ├── fit_none_0.5_0_0_closed_f.log
│       │   ├── fit_sample-dists_0.5_1_1_closed_t.log
│       │   └── segfault.dag.svg
│       └── mice.sensspec.labmeeting
├── logfiles                             # Logfiles generated by Flux & Mothur
│   └── pbs
├── pbs-torque                           # Snakemake profile for automatically submitting PBS jobs on Flux
│   ├── config.yaml                      # Be sure to edit this before submitting jobs
│   ├── pbs-jobscript.sh
│   ├── pbs-status.py
│   └── pbs-submit.py
└── results
    ├── dataset-as-reference             # Results from using the dataset as its own reference
    │   ├── human
    │   ├── marine
    │   ├── mice
    │   ├── soil
    │   └── soil_1000
    │       └── figures
    │           ├── aggregate.sensspec.mcc.full.png
    │           ├── aggregate.sensspec.mcc.iters.png
    │           └── aggregate.sensspec.mcc.png
    ├── silva-as-reference               # Results from using Silva as the reference
    └── workflows                        # Visuals created by Snakemake --dag
        ├── optifit-pipeline.FULL.dag.svg
        ├── optifit-pipeline.config_test.dag.svg
        ├── optifit-pipeline.test.dag.svg
```
