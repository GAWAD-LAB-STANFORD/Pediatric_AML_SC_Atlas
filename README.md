# Pediatric AML single-cell atlas — analysis code

Analysis code and plot-level source data for the pediatric acute myeloid leukemia (AML)
single-cell RNA-seq atlas and surface-immunotherapy-target study. The study defines
prognosis-associated leukemic populations across cytogenetic subtypes and screens their
surface expression against normal tissues to nominate **CD96** (with **ITGAX** and
**TNFRSF4**) as candidate marrow-sparing immunotherapy targets.

> Pang Y, Lyle AG, Cai B, Sanders LM, Gonzalez-Pena V, Aragon AS, Vaske OM, Gawad C.
> *Single-cell profiling and bulk deconvolution reveal prognosis-associated leukemic
> populations and candidate combinatorial marrow-sparing immunotherapy targets in
> pediatric AML.* (In preparation, 2026.)

## Layout

```
Figure_1 … Figure_7          main-figure panel scripts (+ per-folder source_data/ for Figures 5–7)
Figure_S1 … Figure_S27       supplementary-figure panel scripts (+ per-folder source_data/)
source_data/                 the manuscript Source Data, one folder per figure — the canonical
                             plot-level values behind every main and supplementary panel
shared_CD96_pipeline/        shared environment, utilities, and the surfaceome pipeline
README_code_locations.md     narrative map of where each analysis lives
FIGURE_SCRIPT_MAP.md         figure/panel → script → source-data crosswalk
```

Start with **[`FIGURE_SCRIPT_MAP.md`](FIGURE_SCRIPT_MAP.md)** and
**[`README_code_locations.md`](README_code_locations.md)** to find the script behind any
panel, and **[`source_data/README_source_data_locations.md`](source_data/README_source_data_locations.md)**
for the data behind it.

## Where each figure reads its data

- **Figures 5–7 and S1–S27** read the copies in their own folder's `source_data/`.
- **Figures 1–4 and S9** (and shared inputs such as the reference leukemic-state
  assignments in `source_data/Figure_3/v3/`) read from the top-level `source_data/` tree.

Scripts locate the package through a root path: `Figure_3/config.py` defines `BASE` (with a
`FIG3_DATA` environment-variable override for the Figure 3 data directory), and the R scripts
carry the root path near the top of each file. Point these at your clone, with `source_data/`
standing in for the package's `08_Source_Data/` folder. Pinned environments live under
`shared_CD96_pipeline/` (Python and R package snapshots); scripts are a mix of Python (scanpy)
and R.

The raw single-cell atlas itself (the `.h5ad`) is **not** stored here — it is large and lives
in the public data repositories below; everything in `source_data/` is the plot-level data
derived from it.

## Data availability

- **Single-cell gene-expression matrices** — Alex's Lemonade Stand Foundation
  Single-cell Pediatric Cancer Atlas, project **SCPCP000007**.
- **TARGET pediatric AML bulk RNA-seq / outcomes** — dbGaP **phs000465**.
- **Beat AML** adult validation cohort — as published by the Beat AML consortium.
- **Normal-tissue expression** — CZ CELLxGENE Discover Census.

## License

Code is released under the [MIT License](LICENSE). Third-party data accessed through the
repositories above remain under their respective terms.
