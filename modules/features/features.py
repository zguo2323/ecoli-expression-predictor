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
        return 0.0
    tail = promoter_seq[-30:] if len(promoter_seq) > 30 else promoter_seq
    junction = (tail + rbs_seq).upper().replace("T", "U")
    _, mfe = RNA.fold(junction)
    return float(mfe)


def extract_all_features(promoter_seq: str, rbs_seq: str) -> dict:
    spacer = get_spacer_length(promoter_seq)
    sd_spacing = get_sd_spacing(rbs_seq)
    return {
        "gc_promoter": compute_gc_content(promoter_seq),
        "gc_rbs": compute_gc_content(rbs_seq),
        "score_minus10": score_minus10_box(promoter_seq),
        "score_minus35": score_minus35_box(promoter_seq),
        "spacer_length": spacer,
        "spacer_optimal": 15 <= spacer <= 21,
        "score_sd": score_sd_sequence(rbs_seq),
        "sd_spacing": sd_spacing,
        "sd_spacing_optimal": 5 <= sd_spacing <= 10,
        "mrna_folding_energy": compute_mrna_folding_energy(promoter_seq, rbs_seq),
    }
