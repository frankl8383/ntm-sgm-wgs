from scripts.random_pair_sensitivity import merge_intervals, subtract_covered


def test_merge_intervals_removes_overlap():
    assert merge_intervals([(10, 20), (15, 30), (40, 50)]) == [(10, 30), (40, 50)]


def test_subtract_covered_returns_unseen_segments():
    assert subtract_covered((10, 30), [(0, 15), (20, 25)]) == [(15, 20), (25, 30)]
