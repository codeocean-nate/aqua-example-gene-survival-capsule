# Gene expression survival analysis (Kaplan-Meier + Cox)

Does the expression of a gene in a tumor cohort track with patient overall survival? This capsule
answers that question for one gene, or several at once. Give it a gene symbol (for example
`TP53`) and it splits the cohort into high- and low-expression groups, draws a Kaplan-Meier curve,
runs a log-rank test, and fits a Cox proportional-hazards model on the continuous expression value
adjusted for clinical covariates. It is meant for a biologist or bioinformatician who wants a quick
"is this gene prognostic?" check and expects to **rerun it for a new gene** by changing `--gene`.

## Method

1. **Gene lookup.** `--gene` is matched case-insensitively against the `symbol` and `gene_id`
   columns of `synthetic-gene-annotation/genes.tsv` (so `tp53`, `TP53` and `SYN000001` all work).
   The gene must also be one of the 100 genes measured in the cohort; if it is not, the script stops
   before running anything and prints the closest matching cohort genes and the list of named
   genes available.
2. **Cohort assembly.** Expression columns are matched to `clinical.tsv` by `patient_id`; patients
   with a missing value in the time, event or covariate columns are dropped.
3. **High/low groups.** `--split median` puts patients above the gene's median in "High" and the rest
   in "Low". `--split quartile` compares the top quartile (>= 75th percentile) with the bottom
   quartile (<= 25th percentile) and drops the middle 50%.
4. **Kaplan-Meier and log-rank.** `survival::survfit(Surv(time, event) ~ group)` for the curves and
   95% confidence bands; `survival::survdiff` for the log-rank test (chi-square, 1 df) between the
   High and Low groups.
5. **Cox model.** `survival::coxph(Surv(time, event) ~ expression + <covariates>)` on **all**
   patients, with expression as a continuous variable (log2 units). The reported hazard ratio (HR) is
   per 1 log2 unit increase in expression, with a Wald 95% CI and p-value. HR > 1 means higher
   expression goes with worse survival. Character covariates (`stage`, `sex`) are treated as
   categorical factors (stage reference level is `I`); `age` is numeric. Use `--covariates none` for
   an unadjusted model.
6. **Plots and tables** are made with base R + ggplot2 (no other packages). The significance
   statement in `summary.txt` uses p < 0.05 on the Cox p-value; there is no multiple-testing correction
   across genes.

## Inputs

Two data assets must be attached (Code Ocean mounts them read-only under `/data/<name>/`).

| Data asset (mount name) | Path | Files used | Required columns |
|---|---|---|---|
| `synthetic-tumor-cohort` | `/data/synthetic-tumor-cohort/` | `expression.tsv` | `gene_id`, then one column per patient (`P0001`..`P0240`), log2 TPM-like values (100 genes x 240 patients) |
| | | `clinical.tsv` | `patient_id`, `os_months`, `os_event` (1 = death, 0 = censored), `stage` (I-IV), `age`, `sex` |
| `synthetic-gene-annotation` | `/data/synthetic-gene-annotation/` | `genes.tsv` | `gene_id`, `symbol` (other columns are ignored) |

The data location can be overridden with the `DATA_DIR` environment variable (default `/data`);
the analysis reads `$DATA_DIR/synthetic-tumor-cohort/` and `$DATA_DIR/synthetic-gene-annotation/`.

## Parameters

| Flag | Default | Meaning |
|---|---|---|
| `--gene` | `TP53` | Gene symbol (case-insensitive) or gene_id. Comma-separate several, e.g. `TP53,EGFR`; each gene is analysed and all appear in the combined table. |
| `--split` | `median` | How to make the High/Low groups for the KM curve and log-rank test: `median` or `quartile` (top vs bottom quartile, middle dropped). The Cox model always uses continuous expression. |
| `--covariates` | `age,stage` | Comma-separated clinical columns to adjust the Cox model for (`age`, `stage`, `sex`, or any other `clinical.tsv` column), or `none`. |
| `--time-col` | `os_months` | Clinical column holding follow-up time. |
| `--event-col` | `os_event` | Clinical column holding the event indicator (1 = event, 0 = censored). |

Flags can be written `--gene TP53` or `--gene=TP53`. The parameters used are printed at startup and
saved in `run_summary.json`.

## Outputs

Everything is written to `/results` (override with the `RESULTS_DIR` environment variable).

| File | Contents |
|---|---|
| `km_<SYMBOL>.png` | One Kaplan-Meier plot per gene (e.g. `km_TP53.png`): High vs Low curves with 95% bands, censor marks, log-rank p, and Cox HR/CI/p in the annotation. |
| `cox_results.csv` | One row per gene: `gene` (gene_id), `symbol`, `HR`, `CI_low`, `CI_high`, `p` (Cox Wald p for expression), `logrank_p`, `n_high`, `n_low`. |
| `summary.txt` | Plain-text run summary: parameters, cohort size and event count, and per-gene HR/CI/p, log-rank p, group sizes and a one-line reading of direction. |
| `run_summary.json` | Parameters used, genes analysed, patient/event counts and the list of output files. |

## How to run

With defaults (TP53, median split, adjusted for age and stage):

```bash
code/run
```

Another gene, unadjusted:

```bash
code/run --gene ESR1 --covariates none
```

Several genes, comparing top vs bottom quartile:

```bash
code/run --gene EGFR,MYC --split quartile
```

## Environment

- R 4.2 or newer (developed on R 4.6), base R plus:
  - `survival` (>= 3.5)
  - `ggplot2` (>= 3.4)
- Both are listed in `code/R-packages.txt`. No other packages are needed; arguments are parsed with base R.

Synthetic data notice: all data in this capsule (`synthetic-tumor-cohort`, `synthetic-gene-annotation`) is simulated and is not from real patients.
