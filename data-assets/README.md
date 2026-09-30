# Example data for this capsule (synthetic)

This repository is a Code Ocean capsule, *Gene expression survival analysis (Kaplan-Meier + Cox)*.
It's an example for testing [capsule-reference-skill](https://github.com/codeocean-nate/aqua-capsule-reference-skill).
The capsule reads two data assets, and the files for them are in this folder:

| Data asset (name = mount) | Files | Mounted at |
|---|---|---|
| `synthetic-tumor-cohort` | `expression.tsv`, `clinical.tsv`, `DATASET.md` | `/data/synthetic-tumor-cohort` |
| `synthetic-gene-annotation` | `genes.tsv`, `DATASET.md` | `/data/synthetic-gene-annotation` |

All the data is **synthetic**: 240 made-up patients and 100 genes. `_generators/` holds the
deterministic scripts that produced it (`python3 _generators/<name>.py`, which needs numpy).

## Set it up in a deployment

Paste this into a new Aqua chat:

```text
Create a new capsule by copying the contents of this Git repository: https://github.com/codeocean-nate/aqua-example-gene-survival-capsule

Make it an independent capsule, not linked to the repository. Its data-assets/ folder holds two small synthetic example datasets. Create one data asset from each of the subfolders synthetic-tumor-cohort and synthetic-gene-annotation, using the folder name as both the data asset's name and its mount, with that folder's files at the root of the data asset (leave out _generators and README.md). Attach both data assets to the new capsule at /data/synthetic-tumor-cohort and /data/synthetic-gene-annotation. Then delete the data-assets/ folder from the capsule and commit, so the capsule holds only its code, README and environment. Don't run the survival analysis. Don't ask me any questions. Reply with the capsule's full https URL and each data asset's name, ID and full https URL.
```

Then test the skill:

```text
Turn my capsule <capsule URL from the reply above> into a skill, so that next time I or a colleague can rerun this analysis with a small change instead of starting over. Don't ask me any questions: make reasonable choices and tell me what you did.
```

**If Aqua can't create the data assets,** make them in the web UI instead:
1. Upload each folder's files as a new data asset, with the name and mount shown in the table.
2. Attach both data assets to the capsule at those mounts.
3. Delete `data-assets/` from the capsule.
