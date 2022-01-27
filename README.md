# Overview

![badge](https://github.com/NICHD-BSPC/spatial_clustering/actions/workflows/main.yml/badge.svg)

Contact: ryan.dale@nih.gov

Given the position of variants and protein domains within a cDNA, this tool
will run many permutations and calculate empirical p-values for the geometric
mean of pairwise distances between variants (to assess spatial clustering) as
well as for the number of variants falling within each domain.

# Usage

Requires Python 3.6+ along with the NumPy and pyyaml packages. These can be
installed via pip or conda. Tests, if you'd like to run them, also require the
pytest package to be installed.

Edit the `config.yml` file to reflect your cDNA, variants, and domains of
interest (see that file for details). Then run:

```
python domain_and_clustering_analysis.py config.yml
```

The result will be an annotated report printed to stdout.

The default number of permutations is 1e6, but this can be adjusted with the
`--permutations` argument.
