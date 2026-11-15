import json
from pathlib import Path

DNA = set("ATGC")
RNA = set("AUGC")

# RNA/DNA hybrid NN parameters in 100 mM NaCl, corrected values.
# Source: Sugimoto et al., Nucleic Acids Research 48 (2020), correction 49 (2021).
# Values are ΔG°37 in kcal/mol for each RNA dinucleotide / complementary DNA dinucleotide stack.
NN_DG37 = {
    "AA": -0.7, "AC": -1.5, "AG": -1.3, "AU": -0.4,
    "CA": -1.2, "CC": -1.7, "CG": -1.4, "CU": -0.4,
    "GA": -1.5, "GC": -2.0, "GG": -2.3, "GU": -1.4,
    "UA": -0.5, "UC": -1.4, "UG": -1.6, "UU": 0.2,
}
# Initiation contribution at each duplex end (ΔG°37); values in same study.
INIT_GC = 2.0
INIT_AU = 2.6

def validate_dna(seq, expected_len=None):
    if not seq:
        raise ValueError("표적 DNA 서열을 입력하세요.")
    bad = sorted(set(seq) - DNA)
    if bad:
        raise ValueError(f"표적 DNA에는 A/T/G/C만 사용할 수 있습니다. 잘못된 문자: {', '.join(bad)}")
    if expected_len and len(seq) != expected_len:
        raise ValueError(f"표적 DNA는 {expected_len}nt이어야 합니다. 현재 {len(seq)}nt입니다.")
    return seq

def validate_rna(seq, expected_len=None, label="gRNA"):
    if not seq:
        raise ValueError(f"{label} 서열을 입력하세요.")
    bad = sorted(set(seq) - RNA)
    if bad:
        if "T" in bad:
            raise ValueError(f"{label}는 RNA이므로 T가 아니라 U를 사용해야 합니다.")
        raise ValueError(f"{label}에는 A/U/G/C만 사용할 수 있습니다. 잘못된 문자: {', '.join(bad)}")
    if expected_len and len(seq) != expected_len:
        raise ValueError(f"{label}는 {expected_len}nt이어야 합니다. 현재 {len(seq)}nt입니다.")
    return seq

def compare_sequences(dna, guide):
    validate_dna(dna)
    validate_rna(guide)
    if len(dna) != len(guide):
        raise ValueError("표적 DNA와 gRNA의 길이가 같아야 합니다.")
    positions = [i+1 for i, (d, r) in enumerate(zip(dna, guide)) if d == "T" and r == "U" or d == r]
    mismatches = [i+1 for i, (d, r) in enumerate(zip(dna, guide)) if not (d == r or (d == "T" and r == "U"))]
    line = "".join("|" if (d == r or (d == "T" and r == "U")) else "x" for d, r in zip(dna, guide))
    return {
        "guide": guide,
        "mismatch_count": len(mismatches),
        "mismatch_positions": mismatches,
        "pam_proximal_mismatches": sum(1 for p in mismatches if p >= 13),
        "identity": len(positions)/len(dna),
        "match_line": line,
    }

def hybrid_dg37(rna_guide):
    """ΔG°37 for a fully paired RNA/DNA hybrid using 100 mM NaCl NN parameters.
    Only Watson-Crick matched duplexes are accepted. No mismatch penalty is invented.
    """
    validate_rna(rna_guide)
    if len(rna_guide) < 2:
        raise ValueError("NN 계산에는 최소 2nt가 필요합니다.")
    total = 0.0
    for i in range(len(rna_guide)-1):
        dinuc = rna_guide[i:i+2]
        total += NN_DG37[dinuc]
    # End initiation terms: one for each terminal base pair.
    for base in (rna_guide[0], rna_guide[-1]):
        total += INIT_GC if base in "GC" else INIT_AU
    return total

def load_json_upload(upload):
    try:
        data = json.loads(upload.getvalue().decode("utf-8"))
    except Exception as exc:
        raise ValueError(f"CFD JSON을 읽을 수 없습니다: {exc}")
    if not isinstance(data, dict):
        raise ValueError("CFD JSON의 최상위 구조는 객체여야 합니다.")
    return data

def cfd_score_from_json(guide_rna, off_target_dna, pam, mm_scores, pam_scores):
    """Original Doench-style CFD product for one candidate off-target site.

    guide_rna is the designed RNA guide (5'->3'); off_target_dna is the
    candidate off-target protospacer on the non-target DNA strand (5'->3').
    Matrix keys: r<guide-base>:d<complement-of-off-target-base>,<1-based-position>.
    PAM lookup uses the two bases after the variable N in an NGG-style PAM.
    """
    validate_rna(guide_rna)
    validate_dna(off_target_dna)
    if len(off_target_dna) != len(guide_rna):
        raise ValueError("CFD 계산을 위해 gRNA와 off-target protospacer 길이가 같아야 합니다.")
    pam_key = pam[-2:]
    if pam_key not in pam_scores:
        raise ValueError(f"CFD PAM matrix에 '{pam_key}' 항목이 없습니다.")
    score = float(pam_scores[pam_key])
    dna_complement = {"A":"T", "T":"A", "G":"C", "C":"G"}
    for i, (guide_base, dna_base) in enumerate(zip(guide_rna, off_target_dna), 1):
        if (dna_base == "T" and guide_base == "U") or dna_base == guide_base:
            continue
        key = f"r{guide_base}:d{dna_complement[dna_base]},{i}"
        if key not in mm_scores:
            raise ValueError(f"CFD mismatch matrix에 '{key}' 항목이 없습니다.")
        score *= float(mm_scores[key])
    return score
