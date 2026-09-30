"""Generate the synthetic tumor cohort used by the gene-survival-km capsule.

Deterministic (seed 20240611). Stdlib + numpy only.
Output folder: datasets/synthetic-tumor-cohort/
    expression.tsv  100 genes x 240 patients, log2 TPM-like values, keyed by gene_id
    clinical.tsv    patient_id, os_months, os_event, stage, age, sex
    DATASET.md

Gene IDs follow the shared synthetic-gene-annotation dataset: SYN000001..SYN000060 are
the known cancer genes (TP53, EGFR, ...); 40 further genes are drawn from
SYN000061..SYN002000 (they carry SYNG#### symbols in the annotation).

Planted survival associations (per 1 SD of expression, log hazard ratio):
    TP53  +0.42   (high expression -> worse survival)
    MYC   +0.32   (high expression -> worse survival)
    ESR1  -0.45   (high expression -> better survival)
    CDKN2A +0.15, BCL2 -0.12  (weak effects, near the detection limit)
    all other genes: null (independent of survival, apart from chance)
Stage and age also raise hazard; TP53 expression rises slightly with stage so that
covariate adjustment changes the TP53 hazard ratio a little.
"""
import os

import numpy as np

SEED = 20240611
N_PATIENTS = 240
N_GENES = 100
N_KNOWN = 60
ANNOT_MAX = 2000

KNOWN = """TP53 EGFR KRAS BRCA1 BRCA2 MYC PIK3CA PTEN APC RB1 CDKN2A ERBB2 BRAF NRAS
ALK MET ATM CHEK2 IDH1 IDH2 VHL NF1 SMAD4 CTNNB1 NOTCH1 FGFR1 FGFR2 FGFR3 KIT PDGFRA
ESR1 AR GATA3 FOXA1 CDH1 ARID1A KMT2D STK11 KEAP1 MTOR TSC1 TSC2 MDM2 CCND1 CDK4 CDK6
BCL2 MCL1 JAK2 STAT3 VEGFA HIF1A EZH2 DNMT3A TET2 NPM1 FLT3 RUNX1 SF3B1 CD274""".split()
assert len(KNOWN) == N_KNOWN
SYM2ID = {s: f"SYN{i + 1:06d}" for i, s in enumerate(KNOWN)}

# log hazard per 1 SD of expression
EFFECTS = {"TP53": 0.42, "MYC": 0.32, "ESR1": -0.45, "CDKN2A": 0.15, "BCL2": -0.12}
# per-gene expression SD (log2 units) for the planted genes; others drawn at random
FIXED_SD = {"TP53": 1.0, "MYC": 1.2, "ESR1": 1.4, "EGFR": 1.1}
STAGE_EFFECT = 0.35   # log hazard per stage step (I=0 .. IV=3)
AGE_EFFECT = 0.015    # log hazard per year above 62

rng = np.random.default_rng(SEED)

# ---- genes ---------------------------------------------------------------------------
extra_ids = sorted(rng.choice(np.arange(N_KNOWN + 1, ANNOT_MAX + 1), size=N_GENES - N_KNOWN, replace=False))
gene_ids = [f"SYN{i:06d}" for i in range(1, N_KNOWN + 1)] + [f"SYN{i:06d}" for i in extra_ids]
gene_syms = KNOWN + [f"SYNG{i:04d}" for i in extra_ids]

# ---- patients / clinical -------------------------------------------------------------
patients = [f"P{i:04d}" for i in range(1, N_PATIENTS + 1)]
stage_idx = rng.choice(4, size=N_PATIENTS, p=[0.25, 0.30, 0.28, 0.17])
stage = np.array(["I", "II", "III", "IV"])[stage_idx]
age = np.clip(np.round(rng.normal(62, 11, N_PATIENTS)), 30, 90).astype(int)
sex = rng.choice(["F", "M"], size=N_PATIENTS)

# ---- expression (standardised latent z, then scaled to log2 TPM-like values) ---------
Z = rng.standard_normal((N_GENES, N_PATIENTS))
i_tp53, i_myc = gene_syms.index("TP53"), gene_syms.index("MYC")
Z[i_myc] = 0.25 * Z[i_tp53] + np.sqrt(1 - 0.25**2) * Z[i_myc]          # mild TP53-MYC correlation
Z[i_tp53] = Z[i_tp53] + 0.30 * (stage_idx - stage_idx.mean())          # TP53 rises a little with stage
Z[i_tp53] = (Z[i_tp53] - Z[i_tp53].mean()) / Z[i_tp53].std()           # back to unit variance

mu = rng.uniform(3.0, 8.0, N_GENES)
sd = rng.uniform(0.8, 1.4, N_GENES)
for sym, s in FIXED_SD.items():
    sd[gene_syms.index(sym)] = s
expr = np.clip(mu[:, None] + sd[:, None] * Z, 0.0, None)              # log2(TPM+1)-like, floored at 0

# ---- survival ------------------------------------------------------------------------
lp = np.zeros(N_PATIENTS)
for sym, beta in EFFECTS.items():
    lp += beta * Z[gene_syms.index(sym)]
lp += STAGE_EFFECT * stage_idx + AGE_EFFECT * (age - 62)
lp -= 0.35 * 1.6  # roughly centre the stage contribution so the baseline median stays interpretable

shape, median0 = 1.3, 48.0                        # Weibull shape, baseline median (months)
lam = np.log(2) / median0**shape
u = rng.uniform(size=N_PATIENTS)
t_event = (-np.log(u) / (lam * np.exp(lp))) ** (1 / shape)
t_cens = np.minimum(rng.uniform(18, 108, N_PATIENTS), 96.0)   # staggered entry + administrative cut-off
os_months = np.maximum(np.round(np.minimum(t_event, t_cens), 2), 0.1)
os_event = (t_event <= t_cens).astype(int)

# ---- write ---------------------------------------------------------------------------
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "synthetic-tumor-cohort")
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, "expression.tsv"), "w") as f:
    f.write("gene_id\t" + "\t".join(patients) + "\n")
    for g, row in zip(gene_ids, expr):
        f.write(g + "\t" + "\t".join(f"{v:.3f}" for v in row) + "\n")
with open(os.path.join(out, "clinical.tsv"), "w") as f:
    f.write("patient_id\tos_months\tos_event\tstage\tage\tsex\n")
    for i, p in enumerate(patients):
        f.write(f"{p}\t{os_months[i]:.2f}\t{os_event[i]}\t{stage[i]}\t{age[i]}\t{sex[i]}\n")
with open(os.path.join(out, "DATASET.md"), "w") as f:
    f.write(f"""# synthetic-tumor-cohort

**Synthetic data - not real patients or real expression.** Generated by
`_generators/synthetic-tumor-cohort.py` (seed {SEED}).

A small simulated tumor cohort: {N_PATIENTS} patients, {N_GENES} genes. Mount as
`/data/synthetic-tumor-cohort`.

## expression.tsv

Tab-separated, one row per gene, one column per patient.

| column | meaning |
|---|---|
| gene_id | synthetic stable ID (`SYN######`), the key used in `synthetic-gene-annotation/genes.tsv` |
| P0001 ... P{N_PATIENTS:04d} | log2 TPM-like expression (log2(TPM+1)-style, floored at 0) for that patient |

The {N_GENES} genes are SYN000001-SYN000060 (the known cancer genes such as TP53, EGFR, MYC, ESR1)
plus {N_GENES - N_KNOWN} other genes drawn from SYN000061-SYN002000. Gene symbols are **not** in this file; look
them up in `synthetic-gene-annotation`.

## clinical.tsv

Tab-separated, one row per patient.

| column | meaning |
|---|---|
| patient_id | `P0001`..`P{N_PATIENTS:04d}`, matches the expression columns |
| os_months | overall-survival / follow-up time in months |
| os_event | 1 = death observed, 0 = censored |
| stage | tumor stage `I`, `II`, `III`, `IV` |
| age | age at diagnosis, years |
| sex | `F` or `M` |

## Planted structure (ground truth for testing)

High TP53 and high MYC expression are associated with worse survival, high ESR1 with better survival
(CDKN2A and BCL2 have weak effects); all other genes are null apart from chance. Stage and age also
affect survival.
""")
print(f"ok: {expr.shape[0]} genes x {expr.shape[1]} patients, events={os_event.sum()}/{N_PATIENTS}, median follow-up={np.median(os_months):.1f}")
