#!/usr/bin/env python

"""
Permutation testing for spatial clustering of variants and for presence in domains.

Requires NumPy and Python 3.6+
"""

__authors__ = "Arjun Mittal and Ryan Dale"


import argparse
from textwrap import dedent
import numpy as np
import yaml


class SpatialTest(object):
    def __init__(self, variants, domains, length):
        self.variants = sorted(variants)
        self.domains = {
            k: v for k, v in sorted(domains.items(), key=lambda item: item[1])
        }
        self.length = length

    def dist(self, M, normalize=False):
        """
        Return geometric mean of pairwise differences between values of M.

        Parameters
        ----------

        M : list-like
            Integer coordinates representing variants

        cDNA_length : int
            Total length of cDNA, required if normalize=True

        normalize : bool
            Normalize distances by total gene cDNA length, useful if comparing
            between genes.
        """
        M_d = []
        for i in M:
            for j in M:
                if i < j:
                    M_d.append(j - i)

        M_d = np.array(M_d).astype(int)

        if normalize:
            if self.length is None:
                raise ValueError("need to provide cDNA length if normalizing")
            M_d = M_d / self.length

        k = len(M)
        d_k = M_d.prod() ** (1.0 / k)

        # return geometric mean distance
        return d_k

    def permutations(self, n_perms, mutation_number, cDNA_length):
        """
        Returns a list n_perms arrays, each containing a permutation of
        `mutation_number` random variants.

        Parameters
        ----------

        n_perms : int
            Number of permutations to run

        mutation_number : int
            Number of variants to generate on each permutation

        cDNA_length : int
            Length of cDNA; random variants will be selected from a vector of this
            length
        """
        return [
            np.random.RandomState().choice(
                range(1, self.length), size=mutation_number, replace=True
            )
            for x in range(0, n_perms)
        ]

    def domain_location(self, variants, sorted_locations):
        """
        Returns a dictionary with the same keys as `sorted_locations`, but values
        are the number of variants in each domain.

        Requires sorted input.

        Parameters
        ----------

        variants : list of int

        sorted_locations : dict of int
        """
        domain = {key: 0 for key in sorted_locations.keys()}
        for i in variants:
            for j in sorted_locations:
                if i > sorted_locations[j]:
                    continue
                if i <= sorted_locations[j]:
                    domain[j] += 1
                    break
        return domain

    def calculate_dist_p(self, dist_list, actual_dist, n_perms):
        """
        Return the empirical p-value of `actual_dist` found within the list of
        permuted values `dist_list` which is of length `n_perms`
        """
        return (np.sum([actual_dist > dist_perm for dist_perm in dist_list]) + 1) / (
            n_perms + 1
        )

    def calculate_domain_p(
        self, real_location, perm_locations, sorted_domain_locations, n_perms
    ):
        """
        Return a dict, keyed by domains, containing empirical p-values of number of
        observed variants per domain compared to permutations.
        """
        # TODO: only used sorted_domain_locations for keys?
        perm_location_counts = {key: [] for key in sorted_domain_locations.keys()}
        for domain_count_pair in perm_locations:
            for domain in domain_count_pair:
                perm_location_counts[domain].append(domain_count_pair[domain])
        p_dict = {key: 0 for key in sorted_domain_locations.keys()}

        for domain in p_dict:
            p_dict[domain] = (
                np.sum(real_location[domain] < np.array(perm_location_counts[domain]))
                + 1
            ) / (n_perms + 1)
        return p_dict

    def run(self, n_perms):
        """
        Runs permutation testing, using `n_perms` permutations.
        """
        M = np.array(self.variants)
        n = len(M)
        actual_dist = self.dist(M, normalize=True)
        actual_location_count = self.domain_location(M, self.domains)

        mutation_permutations = self.permutations(n_perms, n, self.length)
        perm_dists = [self.dist(perm, normalize=True) for perm in mutation_permutations]
        perm_location_counts = [
            self.domain_location(perm, self.domains) for perm in mutation_permutations
        ]
        dist_p_val = self.calculate_dist_p(perm_dists, actual_dist, n_perms)
        domain_p_val = self.calculate_domain_p(
            actual_location_count, perm_location_counts, self.domains, n_perms
        )

        return dict(
            actual_dist=actual_dist,
            actual_count=actual_location_count,
            dist_p_val=dist_p_val,
            domain_p_val=domain_p_val,
            n_perms=n_perms,
        )

    def report(self, d):
        """
        Print out the results of run()
        """
        contents = dedent(
            f"""\
        length: {self.length}
        permutations: {d['n_perms']}
        actual geometric mean pairwise distance: {d['actual_dist']}
        distance p-value: {d['dist_p_val']}
        number of variants per domain:
        {d['actual_count']}
        p-values per domain:
        {d['domain_p_val']}
        """
        )
        print(contents)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("config", help="Config file")
    ap.add_argument(
        "--permutations",
        "-p",
        help="Number of permutations, default is %(default)s",
        type=int,
        default=1000000,
    )
    args = ap.parse_args()
    config = yaml.load(open(args.config), Loader=yaml.SafeLoader)

    S = SpatialTest(**config)
    results = S.run(n_perms=args.permutations)
    S.report(results)
