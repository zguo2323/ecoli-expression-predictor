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
    build_transcript_prefix_to_start,
    build_transcript_context,
    compute_mrna_folding_energy,
    compute_translation_accessibility,
    load_reporter_cds,
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


def test_build_transcript_prefix_to_start_uses_tss_offset_and_rna_alphabet():
    prefix = build_transcript_prefix_to_start(
        '"AAA CCC"', '"GGG CCA CATATG"', -2
    )
    assert prefix == "CCGGGCCACAUAUG"


def test_build_transcript_prefix_to_start_zero_offset_starts_at_junction():
    prefix = build_transcript_prefix_to_start("AAACCC", "GGGCCACATATG", 0)
    assert prefix == "GGGCCACAUAUG"


@pytest.mark.parametrize(
    "promoter_seq,rbs_seq,tss_best",
    [
        ("", "GGGCCCATATG", 0),
        ("AAACCC", "GGGCCCATATG", 1.5),
        ("AAACCC", "GGGCCCATATG", -7),
        ("AAACCC", "GGGCCCATATAA", 0),
        ("AAACCN", "GGGCCCATATG!", 0),
    ],
)
def test_build_transcript_prefix_to_start_rejects_invalid_inputs(
    promoter_seq, rbs_seq, tss_best
):
    with pytest.raises(ValueError):
        build_transcript_prefix_to_start(promoter_seq, rbs_seq, tss_best)


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


def test_folding_feature_fails_explicitly_without_viennarna(monkeypatch):
    import modules.features.features as feature_module

    monkeypatch.setattr(feature_module, "_VIENNA_AVAILABLE", False)
    with pytest.raises(RuntimeError, match="ViennaRNA is required"):
        feature_module.compute_mrna_folding_energy("ACGT", "ACGT")
    with pytest.raises(RuntimeError, match="ViennaRNA is required"):
        feature_module.compute_translation_accessibility("ACGU", (0, 1), (1, 2))


def test_extract_all_features_keys():
    result = extract_all_features("TTGACATATAATCCGG", "AAAGAGGAGAAA")
    expected_keys = {
        "gc_promoter", "gc_rbs", "score_minus10", "score_minus35",
        "spacer_length", "spacer_optimal", "score_sd",
        "sd_spacing", "sd_spacing_optimal", "sd_unpaired_probability",
        "start_unpaired_probability", "sd_start_opening_energy",
    }
    assert expected_keys == set(result.keys())


def test_extract_all_features_types():
    result = extract_all_features("TTGACATATAATCCGG", "AAAGAGGAGAAA")
    assert isinstance(result["gc_promoter"], float)
    assert isinstance(result["spacer_optimal"], bool)
    assert all(result[name] != result[name] for name in (
        "sd_unpaired_probability", "start_unpaired_probability", "sd_start_opening_energy"
    ))


def test_build_transcript_context_appends_cds_and_coordinates():
    context, sd_interval, start_interval = build_transcript_context(
        "AAACCC", "AAGGAGGAAACATATG", -2, "ATGAAACCCGGG"
    )
    assert context == "CC" + "AAGGAGGAAACAUAUG" + "AAACCCGGG"
    assert context[sd_interval[0]:sd_interval[1]] == "AGGAGG"
    assert context[start_interval[0]:start_interval[1]] == "AUG"


def test_reporter_reference_has_expected_coding_context():
    cds = load_reporter_cds()
    assert len(cds) == 90
    assert cds.startswith("ATG")
    assert all(cds[i:i + 3] not in {"TAA", "TAG", "TGA"} for i in range(0, 90, 3))


def test_translation_accessibility_is_bounded_and_opening_cost_nonnegative():
    values = compute_translation_accessibility(
        "GGGAAAGGAGGAAACCCCAUGGCGCGCGCGCGCGCGC",
        (6, 12), (18, 21),
    )
    assert 0 <= values["sd_unpaired_probability"] <= 1
    assert 0 <= values["start_unpaired_probability"] <= 1
    assert values["sd_start_opening_energy"] >= 0


def test_opening_energy_matches_partition_function_ratio():
    from modules.features.features import build_transcript_context

    context, sd, start = build_transcript_context(
        "TTGACATATAATCCGG", "AAAGAGGAGAAACATATG", -3
    )
    values = compute_translation_accessibility(context, sd, start)
    # This sequence's unconstrained/constrained PFs give approximately 1.44 kcal/mol.
    assert 1.3 < values["sd_start_opening_energy"] < 1.6


def test_extract_all_features_with_tss_computes_accessibility():
    result = extract_all_features("AAACCC", "AAGGAGGAAACATATG", tss_best=-2)
    assert all(result[name] == result[name] for name in (
        "sd_unpaired_probability", "start_unpaired_probability", "sd_start_opening_energy"
    ))
