# Unit tests; run with pytest -v tests.py

from domain_and_clustering_analysis import SpatialTest


def test_domain_count():
    S = SpatialTest(
        variants=[1, 3, 9, 11, 14, 17, 29, 30, 31, 33, 38, 40, 42, 45, 50, 55],
        length=75,
        domains={
            "reg1": 10,
            "reg2": 14,
            "reg3": 19,
            "reg4": 28,
            "reg5": 39,
            "reg6": 50,
        },
    )
    domain_count = S.domain_location(S.variants, S.domains)

    # Expected results (by hand)
    actual_count = {"reg1": 3, "reg2": 2, "reg3": 1, "reg4": 0, "reg5": 5, "reg6": 4}
    assert domain_count == actual_count


def test_distance():
    S = SpatialTest(variants=[1, 8, 10], length=10, domains={"reg1": 10})
    assert S.dist(S.variants) == 5.0132979349645845

    S = SpatialTest(variants=[1, 3], length=10, domains={"reg1": 10})
    assert S.dist(S.variants) == 2 ** 0.5

    # Ensure total length does not affect distances
    S = SpatialTest(variants=[1, 3], length=20, domains={"reg1": 10})
    assert S.dist(S.variants) == 2 ** 0.5

    # Total length should affect pvals
    S1 = SpatialTest(variants=[1, 3], length=10, domains={"reg1": 10}).run(100)
    S2 = SpatialTest(variants=[1, 3], length=20, domains={"reg1": 10}).run(100)
    assert S1["dist_p_val"] != S2["dist_p_val"]

    # But should not affect domain count pvals
    assert S1["domain_p_val"] == S2["domain_p_val"]
