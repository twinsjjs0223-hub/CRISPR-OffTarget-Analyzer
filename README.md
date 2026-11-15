# CRISPR Off-target Analyzer v3

Streamlit-based sequence comparison tool for SpCas9 educational/research support.

## Run locally

```bash
python -m venv .venv
# Windows
.venv\Scripts\activate
# macOS/Linux
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

## Inputs and strand convention

- Target DNA: 20-nt protospacer on the PAM-adjacent **non-target strand**, written 5′→3′.
- gRNA: 20-nt RNA guide, written 5′→3′, using A/U/G/C only.
- The two input strings are aligned in the same written orientation. A perfect match therefore differs only by DNA T versus RNA U.
- The gRNA physically hybridizes to the opposite target DNA strand, not to the non-target strand entered here.
- Optional PAM: DNA triplet immediately following the entered protospacer; canonical SpCas9 PAM is NGG.

## Calculations

### 1. Mismatch profile
Direct position-by-position comparison, numbered 1–20 from the 5′ end of the entered sequences. The interface separately reports mismatches in positions 13–20 as a descriptive PAM-proximal region; this is not claimed to be a universal seed boundary.

### 2. RNA/DNA hybrid thermodynamics
For perfectly paired RNA/DNA hybrids only, ΔG°37 is calculated by summing the 16 RNA/DNA nearest-neighbor ΔG°37 terms and terminal initiation terms reported for 100 mM NaCl:
- Sugimoto et al. (2020), Nucleic Acids Research, DOI: 10.1093/nar/gkaa572
- Correction (2021), DOI: 10.1093/nar/gkab780

Important: these canonical Watson–Crick NN parameters do not provide mismatch-stack thermodynamics. Therefore this implementation deliberately does **not** invent ΔΔG values for mismatched guides. ΔG°37 is reported only for a perfect guide/target match, and it is not a Cas9 cleavage score.

### 3. CFD
CFD is only computed when the original empirical mismatch and PAM matrices are supplied as JSON. No heuristic is mislabeled as CFD. The key convention expected by the function is:
- mismatch key: `r<guide-base>:d<complement-of-off-target-DNA-base>,<1-based-position>`
- PAM key: the two bases following the variable N in the three-base PAM; for `AGG`, the key is `GG`.

The app compares candidate guides with the supplied target protospacer and, if provided, scores additional off-target DNA sites (23 nt = 20-nt protospacer + 3-nt PAM). It does not search a reference genome and cannot identify genome-wide off-target sites. A complete off-target workflow needs candidate off-target protospacers/PAMs from a genome-wide search, then scores each site.

## Sample input
Target DNA:
`GAGTCCGAGCAGAAGAAGAA`

gRNA candidates:
```text
GAGUCCGAGCAGAAGAAGAA
GAGUCCGAGCAGAAGAAGAU
GAGUCCGAGCAGAAGCAGAA
```

## Tests
```bash
python -m unittest discover -s tests -v
```

## Scientific limitations
This is a transparent sequence/thermodynamics teaching tool, not a clinically validated predictor. It does not model chromatin accessibility, cell type, Cas9 concentration, DNA repair, chromatin state, or actual cleavage outcomes. Validate nominated off-targets experimentally (e.g., GUIDE-seq or CIRCLE-seq) and use a genome-wide search tool for real genomic off-target discovery.
