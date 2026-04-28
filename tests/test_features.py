import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from modules.features.features import (
    compute_gc_content,
    score_minus10_box,
    score_minus35_box,
    score_sd_sequence,
    get_spacer_length,
    get_sd_spacing,
    compute_mrna_folding_energy,
    extract_all_features,
)


def test_compute_gc_content_known():
    assert compute_gc_content("GCGC") == 1.0
    assert compute_gc_content("ATAT") == 0.0
    assert compute_gc_content("ATGC") == 0.5


def test_compute_gc_content_empty():
    assert compute_gc_content("") == 0.0


def test_compute_gc_content_case_insensitive():
    assert compute_gc_content("gcgc") == compute_gc_content("GCGC")


def test_score_minus10_perfect_consensus():
    assert score_minus10_box("AAATATAATAGGG") == 1.0


def test_score_minus35_perfect_consensus():
    assert score_minus35_box("AAATTGACAAAA") == 1.0


def test_score_minus10_no_match():
    assert score_minus10_box("GGGGGGGGGGGG") < 0.4


def test_score_sd_perfect_consensus():
    assert score_sd_sequence("AAAGGAGGAAAA") == 1.0


def test_score_sd_no_match():
    assert score_sd_sequence("TTTTTTTTTTTT") < 0.4


def test_get_spacer_length_known():
    # TTGACA at pos 0, TATAAT at pos 23 → spacer = 23 - 6 = 17
    seq = "TTGACA" + "A" * 17 + "TATAAT"
    assert get_spacer_length(seq) == 17


def test_get_sd_spacing_known():
    # AGGAGG at pos 0, then 7 nt → spacing = 7
    seq = "AGGAGG" + "A" * 7
    assert get_sd_spacing(seq) == 7


def test_compute_mrna_folding_energy_returns_float():
    energy = compute_mrna_folding_energy("TTGACATATAATCCGG", "AAAGAGGAGAAATTTA")
    assert isinstance(energy, float)


def test_compute_mrna_folding_energy_context_dependent():
    # Same RBS paired with different promoters should give different MFEs
    energy1 = compute_mrna_folding_energy("TTGACATATAATCCGG", "AAAGAGGAGAAA")
    energy2 = compute_mrna_folding_energy("GCGCGCGCGCGCGCGC", "AAAGAGGAGAAA")
    assert energy1 != energy2


def test_extract_all_features_keys():
    result = extract_all_features("TTGACATATAATCCGG", "AAAGAGGAGAAA")
    expected_keys = {
        "gc_promoter", "gc_rbs", "score_minus10", "score_minus35",
        "spacer_length", "spacer_optimal", "score_sd",
        "sd_spacing", "sd_spacing_optimal", "mrna_folding_energy",
    }
    assert expected_keys == set(result.keys())


def test_extract_all_features_types():
    result = extract_all_features("TTGACATATAATCCGG", "AAAGAGGAGAAA")
    assert isinstance(result["gc_promoter"], float)
    assert isinstance(result["spacer_optimal"], bool)
    assert isinstance(result["mrna_folding_energy"], float)
