# Fig 1 (orig Fig 1) — single-cell atlas + sample dendrogram

Legacy figure carried from the original manuscript analysis (`Older_Files/ALSF AML scRNA-seq Paper/ALSF_AML_plot/`).

`ALSF_AML_Dendrogram_plotting.R`, `Figure_1_Dendrogram_plotting.R` build the sample dendrogram from the PCA-coordinate CSVs (source data). Lineage/cell-type UMAPs were produced in the upstream scanpy pipeline (see MANIFEST — raw AnnData not shipped).

## `patient_characteristics.csv` — the cohort table

One row per sequenced sample (28 AML, 2 healthy marrow): study id, diagnosis and
FAB subtype, sex, age at diagnosis, risk group, CNS status, disease status and
sample type, WBC and blast % at diagnosis, MRD course, cytogenetics/mutations,
karyotype, vital status, and relapse after induction.

This replaces the biobank spreadsheet that previously sat here. That file was a
938-row export of the whole tissue bank, of which only these 30 rows belong to
this study, and it carried material that is not de-identified: exact specimen
collection dates, a named pathologist, and verbatim flow-cytometry and pathology
reports on per-patient sheets. It has been removed from this repository and from
its git history, and is not redistributable.

The table here is Safe Harbor de-identified: no dates finer than the year (the
one date embedded in a free-text field was reduced to its year), no names, no
record numbers, and every age is well under 90. Anything beyond these fields
belongs with the controlled-access deposit, not here.
