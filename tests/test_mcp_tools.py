import json
from pathlib import Path

import server


def test_only_user_facing_tools_are_registered():
    assert set(server.registered) == {
        "compute_gc_content",
        "extract_all_features",
        "get_sd_spacing",
        "get_spacer_length",
        "predict_expression",
        "rank_rbs_for_promoter",
        "score_minus10_box",
        "score_minus35_box",
        "score_sd_sequence",
    }


def test_non_finite_feature_values_serialize_as_json_null():
    class StubTool:
        def run(self, **kwargs):
            return {"accessibility": float("nan"), "nested": [float("inf"), 0.5]}

    tool = server._make_tool_fn(StubTool(), [], "stub", "")

    assert json.loads(tool()) == {"accessibility": None, "nested": [None, 0.5]}


def test_context_sensitive_tools_prompt_for_missing_measured_tss():
    root = Path(__file__).resolve().parents[1]
    for relative_path in (
        "modules/model/predict_expression.json",
        "modules/model/rank_rbs_for_promoter.json",
        "modules/features/extract_all_features.json",
    ):
        config = json.loads((root / relative_path).read_text())
        description = config["description"].lower()
        assert "ask once" in description
        assert "do not infer or invent" in description
