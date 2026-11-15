import streamlit as st
import pandas as pd
import numpy as np
import re

st.set_page_config(page_title="CRISPR Off-target Analyzer", page_icon="🧬", layout="wide")
st.title("🧬 CRISPR Off-target Analyzer")
st.caption("SpCas9 · RNA guide (U) · DNA protospacer comparison · educational prototype")

st.warning(
    "탐구용 프로토타입입니다. 입력 서열만으로 실제 유전체의 모든 오프타겟 위치를 찾거나 "
    "세포 내 절단 여부를 확정하지 않습니다. CFD 및 ΔΔG는 아래에서 데이터 출처와 계산 수준을 구분해 표시합니다."
)

st.markdown("""
**서열 방향 규칙 (중요)**  
- 표적 DNA 입력: PAM이 붙어 있는 **비표적(non-target) 가닥의 protospacer 20 nt**, 5′→3′  
- gRNA 입력: **RNA spacer 20 nt**, 5′→3′, 염기는 A/U/G/C  
- 따라서 완전 일치하는 gRNA는 표적 DNA의 T를 U로 바꾼 서열과 같습니다.  
- 표적 DNA의 반대 가닥을 입력하려면 먼저 역상보 서열로 바꿔야 합니다.
""")

def clean_dna(s):
    s = re.sub(r'\s+', '', str(s)).upper()
    if not s or set(s) - set("ACGT"):
        raise ValueError("DNA 서열에는 A, C, G, T만 사용할 수 있습니다.")
    return s

def clean_rna(s):
    s = re.sub(r'\s+', '', str(s)).upper()
    if not s or set(s) - set("ACGU"):
        raise ValueError("gRNA는 RNA이므로 A, U, G, C만 사용할 수 있습니다. T는 U로 바꿔 입력하세요.")
    return s

# Explicit teaching weights only. Not the published CFD matrix.
# position 20 is PAM-proximal in the convention used by Doench CFD.
def teaching_mismatch_weight(rna, dna, pos):
    base = {
        ("G","T"):0.87, ("G","A"):0.67, ("C","T"):0.87,
        ("A","G"):0.75, ("U","C"):0.70, ("C","A"):0.60,
        ("A","C"):0.55, ("U","G"):0.55, ("G","C"):0.50,
        ("C","G"):0.45, ("A","U"):0.60, ("U","A"):0.60
    }.get((rna,dna), 0.50)
    # Transparent positional modifier for demonstration, not empirical CFD.
    position_factor = 1.0 - 0.012 * (pos-1)
    return max(0.05, min(1.0, base * position_factor))

def analyze(target_dna, guide_rna, pam):
    if len(target_dna) != 20 or len(guide_rna) != 20:
        raise ValueError("현재 버전은 protospacer와 gRNA spacer 각각 20 nt를 요구합니다.")
    # RNA guide pairs with target strand; provided DNA is the non-target strand
    # and therefore has the same sequence as guide after T->U conversion.
    target_as_rna = target_dna.replace("T","U")
    mismatches=[]
    cfd_like=1.0
    ddg_proxy=0.0
    for i,(r,d) in enumerate(zip(guide_rna,target_as_rna), start=1):
        if r != d:
            w=teaching_mismatch_weight(r,d,i)
            cfd_like *= w
            # illustrative mismatch destabilization contribution; not a measured ΔG
            destabilization = {"G:T":0.5,"U:G":0.5}.get(f"{r}:{d}",0.8)
            pos_mult = 0.75 + 0.025*i
            ddg_proxy += destabilization*pos_mult
            mismatches.append({
                "위치 (1= PAM-distal, 20=PAM-proximal)":i,
                "gRNA 염기 (RNA)":r,
                "표적 DNA 염기":d.replace("U","T"),
                "미스매치":f"r{r}:d{d.replace('U','T')}",
                "영역":"PAM-proximal" if i>=11 else "PAM-distal",
                "교육용 위치가중치":round(w,3)
            })
    pam_ok = bool(re.fullmatch(r"[ACGT]GG", pam))
    return {
        "미스매치 수":len(mismatches),
        "미스매치 위치":", ".join(str(x["위치 (1= PAM-distal, 20=PAM-proximal)"]) for x in mismatches) or "없음",
        "PAM":pam,
        "PAM 적합성": "NGG" if pam_ok else "NGG 아님",
        "CFD 대체 지표":round(cfd_like,4) if pam_ok else 0.0,
        "ΔΔG 불안정성 대리지표":round(ddg_proxy,3),
        "상세":mismatches
    }

st.subheader("1. 서열 입력")
col1,col2=st.columns(2)
with col1:
    target_input=st.text_input("표적 DNA protospacer (20 nt, 5′→3′)", "GAGTCCGAGCAGAAGAAGAA")
    pam_input=st.text_input("표적 부위 PAM (3 nt)", "AGG")
with col2:
    guides_input=st.text_area("gRNA 후보군 (RNA, 후보마다 한 줄: 이름,서열)", 
"""gRNA-완전일치,GAGUCCGAGCAGAAGAAGAA
gRNA-GU-예시,GAGUCCGAGUAGAAGAAGAA
gRNA-원거리-1개, AAGUCCGAGCAGAAGAAGAA
gRNA-다중미스매치,GAGUACGAGCAGAU GAAGAA""", height=160)

st.caption("※ 위 예시는 시연용 가상 서열입니다. gRNA 서열에는 T가 아닌 U를 사용합니다.")
if st.button("분석 실행", type="primary", use_container_width=True):
    try:
        target=clean_dna(target_input)
        pam=clean_dna(pam_input)
        if len(pam)!=3: raise ValueError("PAM은 3 nt여야 합니다.")
        records=[]
        for line_no,line in enumerate(guides_input.splitlines(), start=1):
            if not line.strip(): continue
            if "," not in line: raise ValueError(f"{line_no}번째 줄은 '이름,서열' 형식이어야 합니다.")
            name,seq=line.split(",",1)
            guide=clean_rna(seq)
            result=analyze(target,guide,pam)
            records.append({"후보":name.strip(),"gRNA (RNA)":guide,**{k:v for k,v in result.items() if k!="상세"}})
        if not records: raise ValueError("gRNA 후보를 한 개 이상 입력하세요.")
        df=pd.DataFrame(records).sort_values("CFD 대체 지표",ascending=False)
        st.subheader("2. 후보별 비교 결과")
        st.dataframe(df,use_container_width=True,hide_index=True)
        st.bar_chart(df.set_index("후보")[["CFD 대체 지표","ΔΔG 불안정성 대리지표"]])
        st.subheader("3. 미스매치 상세")
        for line_no,line in enumerate(guides_input.splitlines(), start=1):
            if not line.strip(): continue
            name,seq=line.split(",",1)
            result=analyze(target,clean_rna(seq),pam)
            with st.expander(f"{name.strip()} — 미스매치 {result['미스매치 수']}개"):
                if result["상세"]:
                    st.dataframe(pd.DataFrame(result["상세"]),use_container_width=True,hide_index=True)
                else:
                    st.success("표적 DNA와 완전 일치합니다.")
        csv=df.to_csv(index=False).encode("utf-8-sig")
        st.download_button("CSV 결과 다운로드",csv,"crispr_offtarget_comparison.csv","text/csv")
    except Exception as e:
        st.error(str(e))

st.subheader("4. 점수 해석과 모델의 한계")
st.markdown("""
- **CFD 대체 지표:** 실제 Doench et al. CFD가 아닙니다. 공개된 CFD 행렬 전체를 내장하지 않은 상태에서 설명용 가중치로 산출한 비교 지표입니다. 논문 CFD 점수로 인용하면 안 됩니다.
- **ΔΔG 불안정성 대리지표:** 실측 kcal/mol 단위의 자유에너지 변화가 아닙니다. 실제 RNA-DNA nearest-neighbor 파라미터를 적용한 열역학 계산도 아니므로 숫자를 ΔΔG(kcal/mol)로 해석하지 않습니다.
- **PAM:** 입력한 PAM이 NGG인지 확인합니다. 이 앱은 전체 유전체를 검색하지 않으며, 입력한 표적 자리와 후보 gRNA 간 서열 차이만 비교합니다.
- **구조적 활성화:** REC3, L1/L2, HNH, RuvC의 구조 변화나 절단 반응속도론을 직접 계산하지 않습니다.
- **결론:** 결과는 후보 간 교육용 상대 비교이며, 절단 확률(%) 또는 생물학적 안전성 판정이 아닙니다.
""")
st.markdown("**정확한 CFD로 확장하려면:** Doench et al. 2016 Supplementary Table 19의 position-specific RNA:DNA mismatch activity 값을 불러와 각 미스매치 가중치를 곱해야 합니다. CFD는 각 미스매치 위치·염기 조합별 실험 유래 가중치의 곱으로 계산됩니다.")
st.caption("참고: Doench et al., Nature Biotechnology (2016), Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9.")
