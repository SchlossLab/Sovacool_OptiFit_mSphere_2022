# OptiFitAnalysis

## Managing software dependencies

I'm using the [conda](https://conda.io/docs/) package manager to manage dependencies for this project. If you don't already have it, I recommend installing the [Miniconda](https://conda.io/miniconda.html) Python 3 distribution. [Here's a link](https://conda.io/docs/_downloads/conda-cheatsheet.pdf) to a helpful cheatsheet for using conda.

Create an environment for the project:
```
conda create --name OptiFitAnalysis --file environment.txt
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
conda list --explicit > environment.txt
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
