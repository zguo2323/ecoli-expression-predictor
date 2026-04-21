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


def test_build_dataset_row_count():
    df = build_dataset(SD01, SD02, SD03)
    assert len(df) > 10000


def test_build_dataset_required_columns():
    df = build_dataset(SD01, SD02, SD03)
    for col in ["promoter_id", "rbs_id", "promo_seq", "rbs_seq", "RNA", "prot", "deltaG"]:
        assert col in df.columns


def test_build_dataset_no_null_sequences():
    df = build_dataset(SD01, SD02, SD03)
    assert df["promo_seq"].isnull().sum() == 0
    assert df["rbs_seq"].isnull().sum() == 0
