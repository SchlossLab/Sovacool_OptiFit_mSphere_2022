# OptiFitAnalysis

## Project progress

See the [Analysis Roadmap](https://github.com/SchlossLab/OptiFitAnalysis/blob/master/AnalysisRoadmap.md).

## Developer Notes

 I have my editor set to convert tabs to spaces with a tab length of 4.
 It's crucial for this to be consistent within Python & Snakemake files.

## Managing software dependencies

I'm using the [conda](https://conda.io/docs/) package manager to manage dependencies for this project.
If you don't already have it, I recommend installing the [Miniconda](https://conda.io/miniconda.html) Python 3 distribution.
[Here's a link](https://conda.io/docs/_downloads/conda-cheatsheet.pdf) to a helpful cheatsheet for using conda.

Create an environment for the project:
```
conda env create --name OptiFitAnalysis --file environment.yaml
```

Activate the environment before running any code:
```
source activate OptiFitAnalysis
```
Be sure to activate the environment from the login node before submitting jobs on the cluster.

Install new packages with:
```
conda install new_package_name
```

Periodically update the dependencies list:
```
conda list --export > environment.export.txt
```

After ending a session, close the active environment with:
```
source deactivate
```

Almost all dependencies are listed in `environment.export.txt`, which will be installed by conda when you create the environment. The exception to this is the `mothur` program.
If you're a member of the Schloss Lab and you're running this analysis on Flux, you can use the mothur binary here: `/nfs/turbo/schloss-lab/bin/mothur-1.42.0/mothur`.
Otherwise, [download the precompiled binary](https://github.com/mothur/mothur/releases) or
[compile from source](https://github.com/mothur/mothur/blob/master/INSTALL.md).
Be sure to use `mothur` version `1.42.0` or higher.

## Snakemake Configuration

The Snakemake workflow relies on a configuration file called `config.yaml`.
It looks like this:

```
mothur_bin: /nfs/turbo/schloss-lab/bin/mothur-1.42.0/mothur
input_dir: data
datasets:
- soil
- mice
- marine
weights:
- none
- ref-abundance
- sample-abundance
- sample-dists
db_version: v132
subsample_test: False
subsample_size: 1000
iterations: 10
replicates: 10
dataset-as-reference: True
silva-as-reference: True
```

- Set `mothur_bin` to the path to your mothur binary if you're not using the default one on Flux.
- Add or remove samples to the `datasets` list as needed.
- Add or remove weight methods to the `weights` list.
    - Only applies to using the dataset as its own reference.
- To run the workflow with just a subset of the input data for debugging purposes:
    - Set `subsample_test` to `True`.
    - Set `subsample_size` to the number of sequences you want to use.
- `iterations` is the number of times the dataset will be randomly subsampled for each fraction.
    - Only applies to using the dataset as its own reference.
- `replicates` is the number of times OptiClust and OptiFit will run on each sample.
- Set `dataset-as-reference` to `True` to run the workflow in `code/analysis/optifit-dataset-as-ref.smk`.
- Set `silva-as-reference` to `True` to run the workflow in `code/analysis/optifit-silva-ref.smk`

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
snakemake -s path/to/snakefile --forceall
```

#### Dry run

Do a dry run to see which jobs snakemake would run without actually running them:
```
snakemake --dryrun -s path/to/snakefile
```
Before committing changes or submitting jobs to the cluster, test your snakefile for syntax errors with a dry run.

#### On the cluster

Edit the cluster config file `cluster.json` with your email for PBS to send job notifications.
On the cluster, create a PBS script with the following command:
```
snakemake --profile pbs-torque
```
Then submit the PBS script to the cluster with `qsub`. See `code/pbs_scripts` for examples.

Profiles for other cluster systems are available in the [snakemake profiles GitHub](https://github.com/snakemake-profiles/doc).

### Visualizing the DAG

Snakemake creates an image representing the directed acyclic graph (DAG) for a workflow with the following command:
```
snakemake --dag -s path/to/workflow.smk | dot -Tsvg > results/workflows/dag.svg
```

Example DAG for `code/data_processing/get-references.smk`:

![get-references.dag.svg](https://github.com/SchlossLab/OptiFitAnalysis/blob/master/results/workflows/get-references.dag.svg)
