# Overview

![badge](https://github.com/NICHD-BSPC/spatial_clustering/actions/workflows/main.yml/badge.svg)

Contact: ryan.dale@nih.gov

Given the position of variants and protein domains within a cDNA, this tool
will run many permutations and calculate empirical p-values for the geometric
mean of pairwise distances between variants (to assess spatial clustering) as
well as for the number of variants falling within each domain.

# Setup

Requires Python 3.6+ along with the NumPy and pyyaml packages (which are
indicated in the `requirements.txt` file. These can be installed via
[pip](https://docs.python.org/3/installing/index.html) or
[conda](https://docs.conda.io/en/latest/).

E.g. with conda,

```
conda create -n spatial-clustering --file requirements.txt
conda activate spatial-clustering
```

Tests, if you'd like to run them, also require the pytest package to be
installed. These are also run automatically via GitHub Actions on every push to
this repo.

# Configuration

Edit the `config.yml` file to reflect your cDNA, variants, and domains of
interest (see that file for details).

# Usage

Run with the following command:

```
python domain_and_clustering_analysis.py config.yml --permutations 1000
```

The default number of permutations is 1e6, but this can be adjusted with the
`--permutations` argument as in the example above.

The result will be an annotated report printed to stdout, e.g.,

```
length: 2922
permutations: 1000
actual geometric mean pairwise distance: 1.0416633395161507e-05
distance p-value: 0.3806193806193806
number of variants per domain:
{'signal': 0, 'proregion': 1, 'linker1': 1, 'PHM': 4, 'linker2': 3, 'PAL': 7, 'linker3': 0, 'TMD': 0, 'CD': 0}
p-values per domain:
{'signal': 0.2647352647352647, 'proregion': 0.01098901098901099, 'linker1': 0.004995004995004995, 'PHM': 0.6153846153846154, 'linker2': 0.19080919080919082, 'PAL': 0.12787212787212787, 'linker3': 0.5224775224775224, 'TMD': 0.32367632367632365, 'CD': 0.7792207792207793}
```
