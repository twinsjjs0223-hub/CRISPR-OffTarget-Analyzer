import streamlit as st
import pandas as pd
from analysis_core import (
    validate_dna, validate_rna, compare_sequences, hybrid_dg37,
    cfd_score_from_json, load_json_upload
)

st.set_page_config(page_title="CRISPR Off-target Analyzer v3", page_icon="🧬", layout="wide")
st.title("CRISPR Off-target Analyzer v3")
st.caption("SpCas9 교육·연구 보조 분석기 | DNA/RNA 검증 · mismatch 분석 · RNA/DNA NN 열역학 · CFD")

with st.sidebar:
    st.header("분석 설정")
    st.subheader("CFD 매트릭스")
    mm_file = st.file_uploader("공식 CFD mismatch matrix (JSON)", type=["json"], key="mm")
    pam_file = st.file_uploader("공식 CFD PAM matrix (JSON)", type=["json"], key="pam")
    st.caption("원본 Doench/Broad 실험 기반 매트릭스가 모두 제공된 경우에만 CFD를 계산합니다.")

left, right = st.columns([1,1])
with left:
    dna = st.text_input("표적 DNA protospacer (비표적 가닥, 5′→3′)", value="GAGTCCGAGCAGAAGAAGAA").upper().strip()
with right:
    guides_text = st.text_area(
        "gRNA 후보 (각 줄에 하나, RNA 5′→3′)",
        value="GAGUCCGAGCAGAAGAAGAA\nGAGUCCGAGCAGAAGAAGAU\nGAGUCCGAGCAGAAGCAGAA",
        height=130
    )

offtargets_text = st.text_area(
    "선택 입력: 잠재적 off-target DNA site (각 줄에 20nt protospacer + 3nt PAM, 5′→3′)",
    value="",
    height=100,
    placeholder="예: GAGTCCGAGCAGAAGAAGAAGG  (20nt protospacer + NGG)"
)

st.info("방향 규칙: 표적 DNA는 PAM 인접 비표적 가닥의 5′→3′ protospacer, gRNA는 이에 같은 방향으로 대응하는 RNA 서열로 입력합니다. 완전 일치 서열은 DNA의 T를 RNA의 U로 바꾼 형태입니다. 실제 혼성체는 gRNA와 반대편 표적 DNA 가닥 사이에 형성됩니다.")

if st.button("분석 실행", type="primary"):
    try:
        validate_dna(dna, expected_len=20)
        guides = [x.strip().upper().replace(" ", "") for x in guides_text.splitlines() if x.strip()]
        if not guides:
            raise ValueError("gRNA 후보를 한 개 이상 입력하세요.")
        for i, seq in enumerate(guides, 1):
            validate_rna(seq, expected_len=20, label=f"gRNA {i}")

        off_sites = []
        for i, line in enumerate(offtargets_text.splitlines(), 1):
            site = line.strip().upper().replace(" ", "")
            if not site:
                continue
            validate_dna(site, expected_len=23)
            off_sites.append((site[:20], site[20:], f"off-target {i}"))

        mm_data = load_json_upload(mm_file) if mm_file else None
        pam_data = load_json_upload(pam_file) if pam_file else None
        cfd_ready = mm_data is not None and pam_data is not None and bool(off_sites)

        rows, details = [], {}
        for idx, guide in enumerate(guides, 1):
            comp = compare_sequences(dna, guide)
            dg = hybrid_dg37(guide) if comp["mismatch_count"] == 0 else None
            best_cfd, best_site = None, None
            if cfd_ready:
                scored = []
                for off_dna, off_pam, site_name in off_sites:
                    score = cfd_score_from_json(guide, off_dna, off_pam, mm_data, pam_data)
                    scored.append((score, site_name, off_dna + off_pam))
                best_cfd, best_name, best_site = max(scored, key=lambda x: x[0])
                best_site = f"{best_name}: {best_site}"
            rows.append({
                "gRNA": f"gRNA {idx}",
                "서열 (5′→3′)": guide,
                "불일치 수": comp["mismatch_count"],
                "불일치 위치 (5′ 기준)": ", ".join(map(str, comp["mismatch_positions"])) or "없음",
                "PAM 근접 불일치 (13–20)": comp["pam_proximal_mismatches"],
                "일치율 (%)": round(comp["identity"]*100, 1),
                "RNA/DNA ΔG°37 (kcal/mol)": round(dg, 3) if dg is not None else "계산 불가",
                "최고 off-target CFD": round(best_cfd, 4) if best_cfd is not None else "off-target·공식 매트릭스 필요",
                "최고 점수 site": best_site or "—"
            })
            details[f"gRNA {idx}"] = comp

        df = pd.DataFrame(rows)
        st.subheader("분석 결과")
        st.dataframe(df, use_container_width=True, hide_index=True)
        st.download_button("결과 CSV 다운로드", df.to_csv(index=False).encode("utf-8-sig"),
                           "crispr_analysis_v3.csv", "text/csv")

        st.subheader("후보별 불일치 비교")
        chart_df = pd.DataFrame({
            "gRNA": [f"gRNA {i}" for i in range(1, len(guides)+1)],
            "불일치 수": [compare_sequences(dna, g)["mismatch_count"] for g in guides]
        }).set_index("gRNA")
        st.bar_chart(chart_df)

        with st.expander("각 후보의 염기별 비교 보기"):
            for name, comp in details.items():
                st.markdown(f"**{name}**")
                st.code(f"DNA  5′-{dna}-3′\nRNA  5′-{comp['guide']}-3′\n     {comp['match_line']}")
                st.write("불일치 위치:", comp["mismatch_positions"] or "없음")

        st.subheader("계산 방법과 해석 한계")
        st.markdown("""
- **염기 검증:** DNA는 A/T/G/C, gRNA는 A/U/G/C만 허용합니다. gRNA의 T는 오류 처리합니다.
- **불일치 분석:** 입력한 비표적 가닥 DNA와 gRNA를 같은 방향으로 정렬해 염기별 차이를 계산합니다.
- **RNA/DNA 열역학:** 완전 상보적 혼성체에 한해 RNA/DNA 최근접 이웃(NN) 매개변수로 ΔG°₃₇를 합산합니다. mismatch가 포함된 duplex에 canonical NN 표를 적용해 ΔΔG를 추정하지 않습니다.
- **CFD:** 실제 CFD는 설계 gRNA와 off-target DNA protospacer의 불일치 종류·위치별 실험 기반 가중치와 PAM 가중치의 곱입니다. 잠재적 off-target site 23nt 서열과 공식 매트릭스가 모두 있어야 계산합니다.
- **범위:** 이 앱은 사용자가 제공한 off-target 후보만 점수화합니다. 참조 유전체 전체를 검색하지 않으므로 genome-wide off-target site를 발굴하거나 세포 내 절단 확률을 예측하지 않습니다.
""")
    except Exception as e:
        st.error(str(e))

st.divider()
st.caption("연구용 예측 모델이 아닌 서열 비교·열역학 교육용 도구입니다. 실제 off-target 검증에는 유전체 기반 탐색과 GUIDE-seq, CIRCLE-seq 등 실험적 검증이 필요합니다.")
