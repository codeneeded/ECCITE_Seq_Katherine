# ECCITE-Seq of HIV-Infected Twins on Antiretroviral Therapy

Single-cell multimodal (scRNA-seq + ADT + TCR/VDJ) analysis of peripheral immune cells from a pair of monozygotic twins **discordant for HIV-1 infection** and on suppressive antiretroviral therapy (ART). The twin design provides a tightly matched genetic and environmental background, isolating the persistent immunological footprint of HIV infection during ART-mediated viral suppression.

> Multimodal profiling (10x Genomics 5' chemistry with ECCITE-seq) captures, in the same cell:
> - the **transcriptome** (Gene Expression),
> - **surface protein epitopes** (Antibody-Derived Tags / ADT panel), and
> - paired **αβ T-cell receptor clonotypes** (VDJ / TCR-seq).

---

## Background

ECCITE-seq (Expanded CRISPR-compatible CITE-seq) extends CITE-seq by simultaneously profiling transcriptome, surface-protein epitopes, and clonotypes/feature barcodes from the same single cell ([Mimitou et al., *Nat. Methods* 2019](https://doi.org/10.1038/s41592-019-0392-0)). Layering a TCR repertoire on top of paired RNA + ADT enables linking T-cell phenotype to clonal identity — an especially powerful design in HIV, where chronic antigen exposure drives clonal expansion and exhaustion that persist even under successful ART.

Studying this in **monozygotic twins discordant for HIV** controls for host genetics and much of the early life environment, so residual differences in immune state are more confidently attributable to HIV/ART itself rather than donor-to-donor variation.

---

## Experimental design

| | |
|---|---|
| **Subjects** | Monozygotic twin pair, one HIV+ (on ART, virally suppressed) and one HIV− *(fill in suppression status / VL / CD4 if you want)* |
| **Sample type** | *(e.g. PBMCs, sorted T cells — fill in)* |
| **Platform** | 10x Genomics Chromium, 5' v2 chemistry |
| **Libraries per sample** | Gene Expression (GEX), Antibody Capture (ADT, with isotype controls), TCR (VDJ) |
| **Cell hashing** | *(HTOs used? if so, which)* |
| **Aligner / counter** | Cell Ranger `multi` *(fill in version)* |

---

## Repository structure

```
ECCITE_Seq_Katherine/
├── Preprocessing/                       # Per-sample preprocessing inputs / intermediates
├── Scripts/                             # Numbered R analysis pipeline (run in order)
│   ├── 1. Preprocessing_DSBNorm.R       # DSB normalization of ADT using empty droplets + isotype controls
│   ├── 2. Seuratv5_QC.R                 # Seurat v5 QC: %mito, nFeature/nCount filtering
│   ├── 3. Seuratv5_QC2_Doublet_Cell_Cycle.R   # scDblFinder doublet removal + cell-cycle scoring
│   ├── 4. Seuratv5_Isotype_Thresholds.R # Isotype-control–based ADT background thresholds
│   └── ...                              # (downstream integration / WNN / TCR scripts)
├── Integration/                         # Cross-modality / cross-sample integration figures
│   ├── RNA_Integration.png
│   ├── seurat_isotype_ADT_Integration.png
│   └── seurat_isotype_WNN_Clusters.png
├── WNN/                                 # Weighted Nearest Neighbor (RNA + ADT) outputs
│   ├── Cluster_Plots/                   # UMAPs / cluster annotations
│   ├── Cluster_Stats/                   # Cell-type composition per condition
│   ├── Average_Expression/              # Per-cluster average expression matrices
│   ├── Feature_Plot/                    # Marker gene / protein feature plots
│   └── Violin_Plot/                     # Marker expression violins
├── VDJ/
│   └── TCR/                             # scRepertoire-based TCR clonotype analyses
├── ECCITE_Seq_Katherine.Rproj
└── .gitignore
```

---

## Pipeline overview

The analysis is implemented end-to-end in **R / Seurat v5**, with the following stages:

1. **ADT denoising — `1. Preprocessing_DSBNorm.R`**
   DSB normalization ([Mulè et al., *Nat. Commun.* 2022](https://doi.org/10.1038/s41467-022-29356-8)) using empty droplets to model ambient antibody background and isotype controls to correct cell-to-cell technical noise. Produces a denoised ADT assay used for all downstream protein-based analyses.

2. **RNA quality control — `2. Seuratv5_QC.R`**
   Filtering on `nFeature_RNA`, `nCount_RNA`, percent mitochondrial reads, and percent ribosomal reads. Per-sample QC plots written for review before thresholds are locked.

3. **Doublet removal & cell-cycle correction — `3. Seuratv5_QC2_Doublet_Cell_Cycle.R`**
   Heterotypic doublet detection (scDblFinder) and Seurat cell-cycle scoring (`CellCycleScoring`) so that cycling state can be regressed out / monitored downstream.

4. **Isotype-control ADT thresholds — `4. Seuratv5_Isotype_Thresholds.R`**
   Per-marker positivity thresholds derived from isotype-control distributions, used to gate ADT-positive populations consistently across samples.

5. **Integration**
   Per-modality integration across the twin samples (RNA, ADT) — see figures in `Integration/`.

6. **Weighted Nearest Neighbor (WNN) multimodal clustering**
   Joint RNA + ADT embedding via Seurat's WNN ([Hao et al., *Cell* 2021](https://doi.org/10.1016/j.cell.2021.04.048)) — clusters, marker discovery, and per-cluster summaries live in `WNN/`.

7. **TCR repertoire — `VDJ/TCR/`**
   Clonotype-level analyses (clonal expansion, diversity, sharing between twins) layered onto the WNN-defined cell types, e.g. with **scRepertoire**.

---

## Dependencies

Core R packages used across the pipeline:

- [`Seurat`](https://satijalab.org/seurat/) ≥ 5.0 — multimodal single-cell framework
- [`dsb`](https://CRAN.R-project.org/package=dsb) — ADT denoising
- [`scDblFinder`](https://bioconductor.org/packages/scDblFinder/) — doublet detection
- [`scRepertoire`](https://www.borch.dev/uploads/screpertoire/) — TCR clonotype analysis
- `SingleCellExperiment`, `SummarizedExperiment` (Bioconductor)
- `ggplot2`, `patchwork`, `dplyr`, `tidyr` — plotting & data wrangling

Cell Ranger is used upstream for read alignment and UMI/feature counting.


## License

*(Add a license — e.g. MIT for code, CC-BY for figures — or state "All rights reserved" if you'd rather decide later.)*
