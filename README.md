# CRISPR Off-target Analyzer (v2)

## RNA 표기
gRNA spacer는 RNA이므로 A/U/G/C로 입력하며 T는 허용하지 않습니다.

## 서열 방향
표적 DNA는 PAM이 인접한 비표적 가닥(non-target strand)의 20-nt protospacer를 5′→3′로 입력합니다. 이 convention에서 완전 일치 gRNA는 표적 DNA의 T를 U로 바꾼 서열입니다. 반대 가닥을 입력하면 역상보 변환이 필요합니다.

## 실행
```bash
pip install -r requirements.txt
streamlit run app.py
```

## GitHub 업로드
```bash
git init
git add .
git commit -m "Build CRISPR off-target analyzer v2"
git branch -M main
git remote add origin https://github.com/USERNAME/REPOSITORY.git
git push -u origin main
```

## Streamlit 배포
Streamlit Community Cloud에서 저장소 연결 후 main file path를 `app.py`로 지정합니다.

## 과학적 한계
현재 CFD 대체 지표와 ΔΔG 불안정성 대리지표는 교육용 휴리스틱입니다. 실제 CFD 행렬 또는 RNA-DNA nearest-neighbor thermodynamic parameters가 아닙니다. 실제 CFD 점수, kcal/mol, 절단 확률로 해석하면 안 됩니다. 정확한 CFD 구현에는 Doench et al. 2016 Supplementary Table 19의 실험 유래 mismatch weights가 필요합니다.
