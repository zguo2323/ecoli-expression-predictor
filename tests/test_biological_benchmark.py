import pandas as pd

from scripts.biological_benchmark import make_group_splits


def test_group_holdouts_have_no_group_leakage_and_double_unseen():
    rows = [
        {"promoter_id": f"p{p}", "rbs_id": f"r{r}"}
        for p in range(5) for r in range(5)
    ]
    df = pd.DataFrame(rows)
    splits = make_group_splits(df, test_size=0.4, random_state=7)

    train, test = splits["promoter_holdout"]
    assert set(df.iloc[train].promoter_id).isdisjoint(df.iloc[test].promoter_id)
    train, test = splits["rbs_holdout"]
    assert set(df.iloc[train].rbs_id).isdisjoint(df.iloc[test].rbs_id)
    train, test = splits["double_unseen"]
    assert set(df.iloc[train].promoter_id).isdisjoint(df.iloc[test].promoter_id)
    assert set(df.iloc[train].rbs_id).isdisjoint(df.iloc[test].rbs_id)
    assert len(test) > 0


def test_group_splits_are_deterministic():
    df = pd.DataFrame([
        {"promoter_id": f"p{p}", "rbs_id": f"r{r}"}
        for p in range(6) for r in range(6)
    ])
    first = make_group_splits(df, random_state=19)
    second = make_group_splits(df, random_state=19)
    for key in first:
        assert all((a == b).all() for a, b in zip(first[key], second[key]))
