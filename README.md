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

The repository is organized **one folder per figure**, each self-contained with the
scripts that build the panels and the plot-level `source_data/` they read:

```
Figure_1 … Figure_7          main-figure panel scripts + source data
Figure_S1 … Figure_S27       supplementary-figure panel scripts + source data
shared_CD96_pipeline/        shared environment, utilities, and the surfaceome pipeline
README_code_locations.md     narrative map of where each analysis lives
FIGURE_SCRIPT_MAP.md         figure/panel → script → source-data crosswalk
```

Start with **[`FIGURE_SCRIPT_MAP.md`](FIGURE_SCRIPT_MAP.md)** and
**[`README_code_locations.md`](README_code_locations.md)** to find the script behind any
panel.

## Reproducing a figure

Each figure folder regenerates its panels from its own `source_data/`, so no external
download is needed to rebuild the figures. Pinned environments live under
`shared_CD96_pipeline/` (Python and R package snapshots). Scripts are a mix of Python
(scanpy) and R.

The raw single-cell atlas itself (the `.h5ad`) is **not** stored here — it is large and
lives in the public data repositories below; the per-figure `source_data/` files are the
plot-level values derived from it.

## Data availability

- **Single-cell gene-expression matrices** — Alex's Lemonade Stand Foundation
  Single-cell Pediatric Cancer Atlas, project **SCPCP000007**.
- **TARGET pediatric AML bulk RNA-seq / outcomes** — dbGaP **phs000465**.
- **Beat AML** adult validation cohort — as published by the Beat AML consortium.
- **Normal-tissue expression** — CZ CELLxGENE Discover Census.

## License

Code is released under the [MIT License](LICENSE). Third-party data accessed through the
repositories above remain under their respective terms.
