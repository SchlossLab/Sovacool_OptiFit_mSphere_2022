# OptiFitAnalysis

## Managing software dependencies

I'm using the [conda](https://conda.io/docs/) package manager to manage dependencies for this project. If you don't already have it, I recommend installing the [Miniconda](https://conda.io/miniconda.html) Python 3 distribution. [Here's a link](https://conda.io/docs/_downloads/conda-cheatsheet.pdf) to a helpful cheatsheet for using conda.

Create an environment for the project:
```
conda env create --name OptiFitAnalysis --file environment.yaml
```

Activate the environment before running any code:
```
source activate OptiFitAnalysis
```

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

Almost all dependencies are listed in `environment.txt`. The exception is  the `mothur` program. Instead of installing with conda, [download the precompiled binary](https://github.com/mothur/mothur/releases) and append the path to the `PATH` in your `.bash_profile`. Then, run:

```
source ~/.bash_profile
```
I'm currently using `mothur` version `1.41.0`.

## Snakemake Workflows

### Dry run

Do a dry run to see which jobs snakemake would run without actually running them.
```
snakemake --dryrun -s path/to/snakefile
```
Before committing changes or submitting jobs to the cluster, test your snakefile for syntax errors with a dry run.

### Visualizing the DAG

Snakemake creates an image representing the directed acyclic graph (DAG) for a workflow with the following command:
```
snakemake --dag -s path/to/workflow.smk | dot -Tsvg > results/figures/workflows/dag.svg
```

Example DAG for `code/data_processing/get-references.smk`:

![get-references.dag.svg](https://github.com/SchlossLab/OptiFitAnalysis/tree/snakemake-analysis/results/figures/workflows/get-references.dag.svg?sanitize=true)
