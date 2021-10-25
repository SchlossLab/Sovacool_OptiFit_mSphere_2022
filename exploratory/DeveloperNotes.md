# Developer Notes

## Project progress

See the [Analysis Roadmap](https://github.com/SchlossLab/OptiFitAnalysis/blob/master/AnalysisRoadmap.md).

## Whitespace in Python

I have my editor set to convert tabs to spaces with a tab length of 4.
In [Atom](https://atom.io) you can use
[:untabify](https://atom.io/packages/tabs-to-spaces) to convert tabs to spaces.
It's crucial for this to be consistent within Python & Snakemake files.
If you get an error while running snakemake like:
```
Unexpected keyword <word> in rule definition (Snakefile, line <line>)
```
It's likely a whitespace issue.

## Managing software dependencies

I use the [conda](https://conda.io/docs/) package manager to manage software
dependencies.
If you don't already have it, I recommend installing the
[Miniconda](https://conda.io/miniconda.html) Python 3 distribution.
If you're experiencing slowness when solving conda environments, follow the
suggestions [here](https://github.com/bioconda/bioconda-recipes/issues/13774)
and also check out [mamba](https://mamba.readthedocs.io/en/latest/)

### Create a conda environment

If you plan to run this workflow on the GreatLakes HPC or another 64-bit Linux
machine, you can get an exact replica of my conda environment with:
```
conda env create --file config/env.export.yaml
```

Otherwise, run:
```
conda env create --name optifit --file config/env.simple.yaml
```

### Activate

Activate the environment before running any code with:
```
conda activate optifit
```

### Update the environment

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
