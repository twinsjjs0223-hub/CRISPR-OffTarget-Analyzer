
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from io import StringIO

st.set_page_config(page_title="CRISPR Off-target Analyzer", page_icon="🧬", layout="wide")

st.title("🧬 CRISPR Off-target Analyzer")
st.caption("gRNA 후보의 서열 불일치, 열역학적 안정성 추정 및 절단 관련 지표를 비교하는 교육·탐구용 시뮬레이터")

st.info(
    "이 앱은 교육·탐구용 프로토타입입니다. 아래 CFD-like 점수와 ΔΔG는 실제 CFD 모델 또는 "
    "실험적으로 검증된 RNA-DNA 최근접 이웃 계산값이 아닌, 경향을 설명하기 위한 휴리스틱 추정치입니다. "
    "실제 오프타겟 안전성·절단 여부 판단이나 임상·실험 설계에 사용하면 안 됩니다."
)

# Illustrative, intentionally transparent mismatch penalties.
# Values are heuristic teaching parameters, not measured thermodynamic parameters.
BASE_PENALTY = {
    "A": {"A": 1.00, "C": 1.15, "G": 0.95, "T": 0.90},
    "C": {"A": 1.10, "C": 1.00, "G": 0.90, "T": 1.05},
    "G": {"A": 0.95, "C": 0.85, "G": 1.00, "T": 0.70},
    "T": {"A": 0.90, "C": 1.00, "G": 0.65, "T": 1.00},
}
TRANSITIONS = {("A","G"),("G","A"),("C","T"),("T","C")}

def clean_seq(seq):
    return "".join(c for c in seq.upper().replace("U", "T") if c in "ACGT")

def mismatch_details(target, guide):
    if len(target) != len(guide):
        raise ValueError("표적 DNA와 gRNA는 같은 길이여야 합니다.")
    out = []
    n = len(target)
    for i, (t, g) in enumerate(zip(target, guide), start=1):
        if t != g:
            # Position 1 is PAM-distal; position n is PAM-proximal in this app's convention.
            pam_prox = i / n
            if (t, g) in TRANSITIONS:
                pair_penalty = 0.78
            else:
                pair_penalty = 1.0
            # PAM-proximal mismatches receive larger heuristic penalty.
            positional = 0.55 + 0.9 * pam_prox
            penalty = pair_penalty * positional
            ddg = BASE_PENALTY[t][g] * positional
            out.append({
                "position": i, "target_base": t, "guide_base": g,
                "mismatch": f"{t}:{g}", "region": "PAM-proximal" if i > n*0.6 else "PAM-distal",
                "heuristic_penalty": penalty, "heuristic_ddG": ddg
            })
    return out

def analyze(target, guide, pam_ok=True):
    mm = mismatch_details(target, guide)
    n = len(target)
    # CFD-like educational score. Not the published CFD model.
    mismatch_load = sum(x["heuristic_penalty"] for x in mm)
    score = 100 * np.exp(-0.34 * mismatch_load)
    score = float(np.clip(score, 0, 100))
    ddg = float(sum(x["heuristic_ddG"] for x in mm))
    # Relative duplex stability proxy: higher penalty -> less stable than perfect match.
    stability = float(100 * np.exp(-0.28 * ddg))
    if not pam_ok:
        score *= 0.05
    # A separate qualitative cleavage-associated index; not a probability.
    cleavage_index = score * (0.65 + 0.35 * stability / 100)
    return {
        "mismatch_count": len(mm),
        "mismatch_positions": ", ".join(str(x["position"]) for x in mm) if mm else "없음",
        "mismatch_types": ", ".join(x["mismatch"] for x in mm) if mm else "완전 일치",
        "PAM": "NGG (입력 확인)" if pam_ok else "미확인/부적합",
        "CFD-like score": round(score, 2),
        "ΔΔG heuristic (kcal/mol)": round(ddg, 2),
        "상대 안정성 지표": round(stability, 2),
        "절단 관련 지표": round(cleavage_index, 2),
        "_details": mm
    }

def example_data():
    target = "GAGTCCGAGCAGAAGAAGA A".replace(" ", "")
    # 20 nt sequences; guide is entered in DNA alphabet (T) or RNA alphabet (U).
    examples = [
        ("gRNA_완전일치", target),
        ("gRNA_G-T_미스매치", target[:9] + ("A" if target[9] != "A" else "C") + target[10:]),
        ("gRNA_PAM원거리_1개", ("A" if target[1] != "A" else "C") + target[1:]),
        ("gRNA_다중미스매치", target[:4] + ("A" if target[4] != "A" else "C") + target[5:13] + ("T" if target[13] != "T" else "G") + target[14:]),
    ]
    return target, examples

with st.sidebar:
    st.header("분석 설정")
    pam_status = st.selectbox("PAM 조건", ["NGG 존재한다고 가정", "PAM 미확인"], index=0)
    st.caption("현재 버전은 20 nt protospacer 비교를 기본으로 합니다. PAM은 별도 확인값으로 처리합니다.")
    st.divider()
    st.markdown("**서열 입력 규칙**")
    st.write("- 표적 DNA: 5′→3′ 방향의 20 nt protospacer")
    st.write("- gRNA: 표적 DNA와 같은 서열 방향으로 입력 (T 또는 U 허용)")
    st.write("- PAM은 protospacer 뒤쪽의 NGG라고 가정")

tab1, tab2 = st.tabs(["직접 입력", "가상 예제 데이터"])

with tab1:
    col1, col2 = st.columns([1, 1])
    with col1:
        target_input = st.text_input("표적 DNA 서열 (20 nt)", value="GAGTCCGAGCAGAAGAAGAA")
    with col2:
        guides_input = st.text_area(
            "gRNA 후보군 (후보마다 한 줄)",
            value="gRNA_A, GAGTCCGAGCAGAAGAAGAA\ngRNA_B, GAGTCCGAGTAGAAGAAGAA\ngRNA_C, AAGTCCGAGCAGAAGAAGAA",
            height=130
        )
    run_manual = st.button("후보 분석 실행", type="primary", use_container_width=True)

with tab2:
    st.write("실제 측정 데이터가 아닌, 미스매치 위치와 종류에 따른 결과 차이를 확인하기 위한 가상 데이터입니다.")
    if st.button("가상 데이터 불러오기", use_container_width=True):
        ex_target, ex_guides = example_data()
        st.session_state["example_target"] = ex_target
        st.session_state["example_guides"] = ex_guides
    if "example_target" in st.session_state:
        st.code(f"표적 DNA: {st.session_state['example_target']}")
        for name, seq in st.session_state["example_guides"]:
            st.code(f"{name}, {seq}")
        if st.button("가상 데이터 분석 실행", type="primary", use_container_width=True):
            st.session_state["run_example"] = True

def run_analysis(target, lines):
    target = clean_seq(target)
    if len(target) != 20:
        st.error(f"표적 DNA는 현재 20 nt여야 합니다. 입력 길이: {len(target)} nt")
        return
    parsed = []
    for idx, line in enumerate(lines, start=1):
        if not line.strip():
            continue
        if "," in line:
            name, seq = line.split(",", 1)
        else:
            name, seq = f"gRNA_{idx}", line
        seq = clean_seq(seq)
        if len(seq) != 20:
            st.warning(f"{name.strip()}: 길이가 20 nt가 아니어서 제외했습니다 ({len(seq)} nt).")
            continue
        result = analyze(target, seq, pam_status == "NGG 존재한다고 가정")
        result["후보"] = name.strip()
        result["gRNA 서열"] = seq
        parsed.append(result)
    if not parsed:
        st.error("분석 가능한 20 nt gRNA 후보가 없습니다.")
        return
    df = pd.DataFrame(parsed).sort_values("절단 관련 지표", ascending=False)
    st.subheader("분석 결과")
    display_cols = ["후보", "gRNA 서열", "mismatch_count", "mismatch_positions", "mismatch_types",
                    "PAM", "CFD-like score", "ΔΔG heuristic (kcal/mol)", "상대 안정성 지표", "절단 관련 지표"]
    st.dataframe(df[display_cols], use_container_width=True, hide_index=True)
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**후보별 CFD-like 점수**")
        chart = df.set_index("후보")[["CFD-like score"]]
        st.bar_chart(chart)
    with c2:
        st.markdown("**후보별 열역학적 불안정성 추정치**")
        chart2 = df.set_index("후보")[["ΔΔG heuristic (kcal/mol)"]]
        st.bar_chart(chart2)
    st.subheader("미스매치 위치 분석")
    for _, row in df.iterrows():
        with st.expander(f"{row['후보']} — 미스매치 {row['mismatch_count']}개"):
            if not row["_details"]:
                st.success("표적 DNA와 완전히 일치합니다.")
            else:
                detail_df = pd.DataFrame(row["_details"])
                st.dataframe(detail_df[["position", "target_base", "guide_base", "mismatch", "region",
                                        "heuristic_penalty", "heuristic_ddG"]], hide_index=True, use_container_width=True)
    csv = df[display_cols].to_csv(index=False).encode("utf-8-sig")
    st.download_button("분석 결과 CSV 다운로드", csv, "crispr_offtarget_results.csv", "text/csv")
    st.caption("해석 주의: ΔΔG heuristic은 교육용 가중치 합산값이며 실측 자유에너지 변화가 아닙니다. CFD-like score 역시 논문 원본 CFD 점수가 아닙니다. 두 지표 모두 절단 확률(%)이 아닙니다.")

if run_manual:
    run_analysis(target_input, guides_input.splitlines())

if st.session_state.get("run_example", False):
    ex_target = st.session_state["example_target"]
    ex_lines = [f"{name}, {seq}" for name, seq in st.session_state["example_guides"]]
    run_analysis(ex_target, ex_lines)
    st.session_state["run_example"] = False

st.divider()
st.markdown("**모델의 범위와 한계**")
st.write(
    "이 프로토타입은 표적 DNA와 gRNA의 서열 정렬을 시연합니다. 실제 GUIDE-seq/Digenome-seq 결과, "
    "검증된 CFD lookup table, RNA-DNA nearest-neighbor thermodynamic parameters, "
    "Cas9 구조·동역학 데이터는 포함하지 않습니다. 따라서 실험적 오프타겟 검증을 대체하지 않습니다."
)
