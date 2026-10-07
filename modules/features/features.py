from functools import lru_cache
from pathlib import Path

try:
    import RNA
    _VIENNA_AVAILABLE = True
except ImportError:
    _VIENNA_AVAILABLE = False


def compute_gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    seq = seq.upper()
    return sum(1 for c in seq if c in "GC") / len(seq)


def _normalize_dna_sequence(seq: str) -> str:
    return "".join(str(seq).replace('"', "").upper().split())


def build_transcript_prefix_to_start(
    promoter_seq: str,
    rbs_seq: str,
    tss_best: int,
) -> str:
    """Build the transcribed sequence from the measured TSS through the start codon.

    ``tss_best`` is the TSS offset from the promoter/RBS junction reported by
    Kosuri et al. The returned sequence intentionally ends at the start codon:
    the reporter CDS must be supplied from a verified construct sequence before
    computing features that depend on downstream coding context.
    """
    promoter = _normalize_dna_sequence(promoter_seq)
    rbs = _normalize_dna_sequence(rbs_seq)

    if not promoter or not rbs:
        raise ValueError("promoter_seq and rbs_seq must be non-empty DNA sequences")
    if any(base not in "ACGTN" for base in promoter + rbs):
        raise ValueError("promoter_seq and rbs_seq may contain only A, C, G, T, or N")
    try:
        tss_value = float(tss_best)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("tss_best must be an integer offset from the promoter/RBS junction") from exc
    if not tss_value.is_integer():
        raise ValueError("tss_best must be an integer offset from the promoter/RBS junction")
    tss_offset = int(tss_value)

    sequence = promoter + rbs
    tss_index = len(promoter) + tss_offset
    start_index = len(sequence) - 3
    if not 0 <= tss_index <= start_index:
        raise ValueError("TSS must fall within the construct and upstream of its terminal start codon")
    if sequence[start_index:] != "ATG":
        raise ValueError("rbs_seq must end with the construct's ATG start codon")

    return sequence[tss_index:start_index + 3].replace("T", "U")


@lru_cache(maxsize=1)
def load_reporter_cds() -> str:
    """Load the provisional 90-nt reporter context documented in data/reference."""
    fasta = Path(__file__).resolve().parents[2] / "data/reference/sfgfp_first90nt.fasta"
    lines = [line.strip() for line in fasta.read_text().splitlines() if line.strip()]
    sequence = _normalize_dna_sequence("".join(line for line in lines if not line.startswith(">")))
    if len(sequence) != 90 or not sequence.startswith("ATG"):
        raise ValueError("Reporter CDS reference must be 90 nt and begin with ATG")
    if any(sequence[i:i + 3] in {"TAA", "TAG", "TGA"} for i in range(0, len(sequence), 3)):
        raise ValueError("Reporter CDS reference contains an in-frame stop codon")
    return sequence


def build_transcript_context(
    promoter_seq: str,
    rbs_seq: str,
    tss_best: int,
    reporter_cds_seq: str | None = None,
) -> tuple[str, tuple[int, int], tuple[int, int]]:
    """Return RNA from measured TSS through 90 nt of CDS and SD/start intervals.

    Intervals are 0-based, half-open coordinates on the returned RNA sequence.
    The reporter CDS begins with the same ATG already at the end of ``rbs_seq``;
    that codon is included only once in the context.
    """
    promoter = _normalize_dna_sequence(promoter_seq)
    rbs = _normalize_dna_sequence(rbs_seq)
    cds = _normalize_dna_sequence(reporter_cds_seq or load_reporter_cds())
    if any(base not in "ACGT" for base in cds):
        raise ValueError("reporter_cds_seq may contain only A, C, G, and T")
    if not cds.startswith("ATG") or len(cds) < 6:
        raise ValueError("reporter_cds_seq must begin with ATG and include downstream coding bases")
    prefix = build_transcript_prefix_to_start(promoter, rbs, tss_best)
    rna = prefix + cds[3:].replace("T", "U")
    tss = int(float(tss_best))
    rbs_start = -tss
    sd_offset = max(range(max(1, len(rbs) - 5)), key=lambda i: sum(a == b for a, b in zip(rbs[i:i + 6], "AGGAGG")))
    sd_start = rbs_start + sd_offset
    start_start = rbs_start + len(rbs) - 3
    return rna, (sd_start, sd_start + 6), (start_start, start_start + 3)


def score_pwm(seq: str, consensus: str) -> float:
    seq = seq.upper()
    consensus = consensus.upper()
    n = len(consensus)
    if len(seq) < n:
        return 0.0
    best = 0.0
    for i in range(len(seq) - n + 1):
        window = seq[i:i + n]
        score = sum(1 for j in range(n) if window[j] == consensus[j]) / n
        if score > best:
            best = score
    return best


def score_minus10_box(promoter_seq: str) -> float:
    return score_pwm(promoter_seq, "TATAAT")


def score_minus35_box(promoter_seq: str) -> float:
    return score_pwm(promoter_seq, "TTGACA")


def get_spacer_length(promoter_seq: str) -> int:
    seq = promoter_seq.upper()
    consensus_35 = "TTGACA"
    consensus_10 = "TATAAT"
    n = 6

    best_35_score, pos_35 = 0.0, -1
    for i in range(len(seq) - n + 1):
        window = seq[i:i + n]
        score = sum(1 for j in range(n) if window[j] == consensus_35[j]) / n
        if score > best_35_score:
            best_35_score = score
            pos_35 = i

    best_10_score, pos_10 = 0.0, -1
    for i in range(len(seq) - n + 1):
        if i <= pos_35:
            continue
        window = seq[i:i + n]
        score = sum(1 for j in range(n) if window[j] == consensus_10[j]) / n
        if score > best_10_score:
            best_10_score = score
            pos_10 = i

    if pos_35 == -1 or pos_10 == -1:
        return -1
    spacer = pos_10 - (pos_35 + n)
    return spacer if spacer >= 0 else -1


def score_sd_sequence(rbs_seq: str) -> float:
    return score_pwm(rbs_seq, "AGGAGG")


def get_sd_spacing(rbs_seq: str) -> int:
    seq = rbs_seq.upper()
    consensus = "AGGAGG"
    n = len(consensus)
    if len(seq) < n:
        return -1

    best_score, pos_sd = 0.0, -1
    for i in range(len(seq) - n + 1):
        window = seq[i:i + n]
        score = sum(1 for j in range(n) if window[j] == consensus[j]) / n
        if score > best_score:
            best_score = score
            pos_sd = i

    if pos_sd == -1:
        return -1
    return len(seq) - (pos_sd + n)


def compute_mrna_folding_energy(promoter_seq: str, rbs_seq: str) -> float:
    if not _VIENNA_AVAILABLE:
        raise RuntimeError("ViennaRNA is required for RNA folding features")
    tail = promoter_seq[-30:] if len(promoter_seq) > 30 else promoter_seq
    junction = (tail + rbs_seq).upper().replace("T", "U")
    _, mfe = RNA.fold(junction)
    return float(mfe)


def compute_translation_accessibility(
    mrna_context: str,
    sd_interval: tuple[int, int],
    start_interval: tuple[int, int],
    temperature: float = 37.0,
    max_bp_span: int = 120,
) -> dict:
    """Compute ensemble unpaired probabilities and SD-through-start opening ΔG."""
    if not _VIENNA_AVAILABLE:
        raise RuntimeError("ViennaRNA is required for translation accessibility features")
    if temperature != 37.0:
        raise ValueError("Translation accessibility is currently defined at fixed 37.0°C")
    sequence = mrna_context.upper().replace("T", "U")
    if not sequence or any(base not in "ACGU" for base in sequence):
        raise ValueError("mrna_context must be a non-empty RNA/DNA sequence using A, C, G, U, or T")
    n = len(sequence)
    for name, interval in (("SD", sd_interval), ("start codon", start_interval)):
        if len(interval) != 2 or not 0 <= interval[0] < interval[1] <= n:
            raise ValueError(f"{name} interval is outside mrna_context")

    opening_start = min(sd_interval[0], start_interval[0])
    opening_end = max(sd_interval[1], start_interval[1])
    md = RNA.md(temperature=temperature, max_bp_span=max_bp_span)
    free = RNA.fold_compound(sequence, md)
    free_pf = free.pf()
    bpp = free.bpp()

    def mean_unpaired(interval: tuple[int, int]) -> float:
        probabilities = []
        for zero_based in range(interval[0], interval[1]):
            pos = zero_based + 1
            paired = sum(
                bpp[min(pos, partner)][max(pos, partner)]
                for partner in range(1, n + 1) if partner != pos
            )
            probabilities.append(max(0.0, min(1.0, 1.0 - paired)))
        return float(sum(probabilities) / len(probabilities))

    constrained = RNA.fold_compound(sequence, md)
    for pos in range(opening_start + 1, opening_end + 1):
        constrained.hc_add_up(pos)
    constrained_pf = constrained.pf()
    opening_energy = max(0.0, float(constrained_pf[1]) - float(free_pf[1]))
    return {
        "sd_unpaired_probability": mean_unpaired(sd_interval),
        "start_unpaired_probability": mean_unpaired(start_interval),
        "sd_start_opening_energy": opening_energy,
    }


class FeaturesMCP:
    def __init__(self, config):
        self.config = config

    def initiate(self):
        pass

    def run(self, **kwargs):
        mcp_name = self.config.get("execution_details", {}).get("mcp_name")
        dispatch = {
            "compute_gc_content":          lambda: compute_gc_content(kwargs["seq"]),
            "score_minus10_box":           lambda: score_minus10_box(kwargs["promoter_seq"]),
            "score_minus35_box":           lambda: score_minus35_box(kwargs["promoter_seq"]),
            "get_spacer_length":           lambda: get_spacer_length(kwargs["promoter_seq"]),
            "score_sd_sequence":           lambda: score_sd_sequence(kwargs["rbs_seq"]),
            "get_sd_spacing":              lambda: get_sd_spacing(kwargs["rbs_seq"]),
            "compute_mrna_folding_energy": lambda: compute_mrna_folding_energy(kwargs["promoter_seq"], kwargs["rbs_seq"]),
            "extract_all_features":        lambda: extract_all_features(kwargs["promoter_seq"], kwargs["rbs_seq"], kwargs.get("tss_best")),
        }
        if mcp_name not in dispatch:
            raise ValueError(f"Unknown mcp_name: {mcp_name}")
        return dispatch[mcp_name]()


def extract_all_features(promoter_seq: str, rbs_seq: str, tss_best: int | None = None) -> dict:
    promoter_seq = _normalize_dna_sequence(promoter_seq)
    rbs_seq = _normalize_dna_sequence(rbs_seq)
    spacer = get_spacer_length(promoter_seq)
    sd_spacing = get_sd_spacing(rbs_seq)
    features = {
        "gc_promoter": compute_gc_content(promoter_seq),
        "gc_rbs": compute_gc_content(rbs_seq),
        "score_minus10": score_minus10_box(promoter_seq),
        "score_minus35": score_minus35_box(promoter_seq),
        "spacer_length": spacer,
        "spacer_optimal": 15 <= spacer <= 21,
        "score_sd": score_sd_sequence(rbs_seq),
        "sd_spacing": sd_spacing,
        "sd_spacing_optimal": 5 <= sd_spacing <= 10,
    }
    if tss_best is None:
        features.update({
            "sd_unpaired_probability": float("nan"),
            "start_unpaired_probability": float("nan"),
            "sd_start_opening_energy": float("nan"),
        })
    else:
        context, sd_interval, start_interval = build_transcript_context(
            promoter_seq, rbs_seq, tss_best
        )
        features.update(compute_translation_accessibility(context, sd_interval, start_interval))
    return features
