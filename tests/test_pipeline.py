import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from modules.pipeline.pipeline import (
    load_promoter_table,
    load_rbs_table,
    load_construct_table,
    build_dataset,
)

SD01 = "data/raw/sd01.xls"
SD02 = "data/raw/sd02.xls"
SD03 = "data/raw/sd03.xls"


@pytest.fixture(scope="module")
def built_dataset(tmp_path_factory):
    output_path = tmp_path_factory.mktemp("transcript-context") / "constructs.parquet"
    df = build_dataset(SD01, SD02, SD03, output_path=str(output_path))
    return df, output_path


def test_load_promoter_table_columns():
    df = load_promoter_table(SD01)
    for col in ["promoter_id", "sequence", "mean_RNA", "mean_prot"]:
        assert col in df.columns


def test_load_promoter_table_no_nulls_in_sequence():
    df = load_promoter_table(SD01)
    assert df["sequence"].isnull().sum() == 0


def test_load_promoter_table_no_quotes():
    df = load_promoter_table(SD01)
    for val in df["promoter_id"]:
        assert not str(val).startswith('"'), f"Quoted value found: {val}"


def test_load_construct_table_filters_bad_prot():
    df = load_construct_table(SD03)
    assert len(df) > 0
    # All bad constructs should be filtered — can't check original bad.prot here
    # but column should not exist in result
    assert "bad.prot" not in df.columns


def test_build_dataset_row_count(built_dataset):
    df, _ = built_dataset
    assert len(df) > 10000


def test_build_dataset_required_columns(built_dataset):
    df, output_path = built_dataset
    for col in [
        "promoter_id", "rbs_id", "promo_seq", "rbs_seq", "RNA", "prot",
        "deltaG", "TSS_best", "TSS_pct_best", "transcript_prefix_to_start",
        "transcript_context",
    ]:
        assert col in df.columns
    assert output_path.exists()


def test_build_dataset_no_null_sequences(built_dataset):
    df, _ = built_dataset
    assert df["promo_seq"].isnull().sum() == 0
    assert df["rbs_seq"].isnull().sum() == 0
    assert df["TSS_best"].notna().all()
    assert df["TSS_pct_best"].notna().all()
    assert df["transcript_prefix_to_start"].notna().all()
    assert df["transcript_prefix_to_start"].str.endswith("AUG").all()
    assert df["transcript_prefix_to_start"].str.fullmatch("[ACGUN]+").all()
    assert df["transcript_context"].notna().all()
    assert df["transcript_context"].str.fullmatch("[ACGUN]+").all()
    assert (df["transcript_context"].str.len() == df["transcript_prefix_to_start"].str.len() + 87).all()
