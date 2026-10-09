# Reproducing the figures

## Layout

Unpack both archives into the **same parent directory**, so they sit side by side:

```
<anywhere>/
  07_Code/         <- 14_Code.zip
  08_Source_Data/  <- 13_Source_Data.zip
  raw/             <- you create this only for the builder scripts (see below)
  ref/             <- small reference files, created only if a builder needs one
  04_Main_Figures/         <- renderers write here
  05_Supplementary_Figures/ <- renderers write here
```

Every path in every script is **relative to this layout**. There are no absolute
paths left in the shipped code.

## Software dependencies

Install these before running the reproduction script. The script checks for them
and names anything missing.

**R** (tested on 4.5.2)

```
install.packages(c("ggplot2","dplyr","scales","patchwork","tidyr","ggrepel",
                   "readr","ggtext","ggnewscale","forcats","png","ggdendro",
                   "flexmix","reshape"))
```

**Python** (tested on 3.14)

```
pip install pandas numpy scipy matplotlib anndata h5py pymupdf seaborn scanpy wordcloud
```

`pymupdf` provides the `fitz` module that the composite scripts import.
`h5py` and `anndata` are needed only by the builder scripts that read the raw
atlas, not by the renderers.

## One command

```
bash 07_Code/reproduce_all_figures.sh            # every figure
bash 07_Code/reproduce_all_figures.sh Figure_5   # one figure
```

It runs two sweeps, because some panels consume PNGs that other panels write,
and prints a summary of what reproduced, what is gated on the raw atlas, and
what failed. Per-script logs land in `07_Code/reproduce_logs/`.

## Two kinds of script

**Renderers** (`fig*.R`, `panel_*`, `composite_*`) rebuild a figure from
`08_Source_Data/` and the per-figure `source_data/` folders that ship inside `07_Code/`.
They need nothing else. Run one as:

```
cd 07_Code/Figure_5/scripts && Rscript fig5.R      # writes ../figures/
```

**Builders** (`build_*`, `compute_*`, `cellchat_*`, `explore_*`) regenerate the
plot-level tables from the raw single-cell atlas. They read from `raw/` and
`ref/`, which are not in these archives because of their size.

## Getting the large inputs (only needed for builders)

| put in | what | where to get it |
|---|---|---|
| `raw/` | `ALSF_AML_Combo_3500_with_PAC_new.h5ad` (45 GB processed atlas) | from the authors on request |
| `raw/` | per-sample matrices | Alex's Lemonade Stand Single-cell Pediatric Cancer Atlas, project **SCPCP000007** |
| `raw/` | `NCBI_TARGET_Star_counts_filtered.csv`, TARGET clinical | NCI TARGET via the Genomic Data Commons, **dbGaP phs000465** |
| `ref/` | `surfaceome_genes.txt` and other small reference lists | SURFY / Cell Surface Protein Atlas, as cited in Methods |

The processed tables behind the interactive browser are also downloadable at
<https://gawadlab.org/downloads/aml-atlas-browser-data.tar.gz>, and the browser
itself runs at <https://gawadlab.org/aml-atlas/>.

Code is also maintained at
<https://github.com/GAWAD-LAB-STANFORD/Pediatric_AML_SC_Atlas>.

## Known limits

- One reference in `shared_CD96_pipeline/build_AllFigures_pdf.py` could not be
  mapped to anything in these archives and is marked
  `SET_THIS_PATH_see_README_REPRODUCIBILITY.md`. It points at a manuscript
  folder, not at data, and no figure depends on it.
- `Figure_1/ALSF_AML_Dendrogram_plotting.R` is listed in the project's own
  `shared_CD96_pipeline/REPRODUCIBILITY_CHECKLIST.md` as containing an
  undefined-variable scratch line. That checklist records the full set of
  outstanding items and ships with this archive.
- Renderers reproduce from the shipped data. End-to-end reproduction from raw
  reads additionally needs the upstream Cell Ranger and scanpy processing
  parameters, which are listed as outstanding in the same checklist.
