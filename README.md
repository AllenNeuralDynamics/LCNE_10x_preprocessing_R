# LCNE 10x snRNA-seq AND retroseq (smarstseq) Preprocessing Pipeline

**Locus Coeruleus Norepinephrine (LC-NE) Cell Identification and Analysis**

This capsule performs sequential preprocessing of 10x Genomics single-nucleus RNA-seq data to identify, filter, and characterize locus coeruleus norepinephrine (LC-NE) neurons using Seurat 5 and Harmony batch correction. The pipeline also includes retro-seq data conversion for cross-modality integration.

**note the entire pipeline runs >1 hours**

## Overview

The pipeline implements a multi-stage workflow:
1. Quality control and Harmony batch correction on full dataset
2. LC-NE population identification via marker gene expression
3. LC-specific subclustering and outlier removal
4. Publication-quality visualization
5. Export to AnnData (h5ad) format for Python workflows
6. Retro-seq data conversion for cross-modality integration

## Input Data

### Downloading the Data

Before running this capsule, you will need to download the required datasets and attach them as data assets, From NEMO!

**snRNAseq data:**
> Download link: *(to be filled)*


**Retro-seq data:**
> Download link: *(to be filled)*

Once downloaded, attach each dataset to this capsule via **Capsule Settings → Data Assets → Add Data Asset**, and ensure the mount names match those listed in the section below.

### Required Data Assets

The pipeline requires two attached data assets:

- **LCNE_10x** (`/data/LCNE_10x/`): Primary 10x snRNA-seq count matrices
- **LCv2** (`/data/LCv2/`): Secondary dataset for integration

Expected file structure:
- Raw count matrices (`.h5`, `.mtx`, or Seurat objects)
- Cell metadata including sample identifiers and QC metrics

## Pipeline Workflow

### Script 1: Build Full Object (`01_build_full_object.R`)

**Purpose:** Initial QC and batch correction of the full dataset

**Steps:**
- Load raw 10x data from both data assets
- Quality control filtering (detected gene counts, mitochondrial %, ribosomal %, doublet scores)
- Counts-per-10k (CP10k) normalization followed by log-normalization
- Harmony batch correction across samples
- Initial dimensionality reduction (PCA, UMAP)
- Preliminary clustering

**Output:** `LC_clustered_harmony.rds` (Seurat object) -> this creates ~ 6GB of data. 

### Script 2: Define LC Population (`02_define_LC_population.R`)
**Purpose:** Identify LC-NE neurons based on marker genes

**Steps:**
- Load batch-corrected Seurat object
- Screen clusters for LC-NE markers (Dbh, Th, Slc6a2, Slc18a2)
- Filter cells expressing canonical LC-NE markers
- Subset to LC-enriched clusters

**Output:** Updated Seurat object with LC population defined

### Script 3: LC Subclustering (`03_LC_subclustering.R`)
**Purpose:** Refine LC-NE population and remove outliers

**Steps:**
- Re-cluster LC-NE subset at higher resolution
- Identify and remove technical outliers/doublets
- Differential expression analysis within LC subtypes
- Final QC and validation

**Output:** Refined LC-NE Seurat object

### Script 4: Publication Figures (`02_b_publication_fig.R`)
**Purpose:** Generate publication-quality visualizations

**Figures produced:**
- UMAP plots colored by cluster, sample, and marker expression
- Violin plots of key LC-NE markers
- Cluster composition and proportion plots
- Quality control metrics

**Output:** Figures saved to `/results/figures/`

### Script 5: Export to h5ad (`04_export_LC_to_h5ad.R`)
**Purpose:** Convert LC subset to AnnData format

**Steps:**
- Extract LC-NE cells from Seurat object
- Convert to h5ad using SeuratDisk or zellkonverter
- Preserve metadata and embeddings

**Output:** `LC_subset.h5ad` for scanpy/Python analysis

### Script 6: Retro-seq Conversion (`05_retroseeq_convert.r`)
**Purpose:** Convert and integrate retro-seq data for cross-modality analysis

**Steps:**
- Load retro-seq input data
- Reformat and align to LC-NE dataset conventions
- Export converted output for downstream integration

**Output:** Converted retro-seq data saved to `/results/`

## Key Parameters

**Quality Control:**
```r
min_genes    <- 2000      # Minimum detected genes per nucleus
max_genes    <- 15000     # Maximum detected genes per nucleus
max_mt_pct   <- 4         # Max mitochondrial % (excludes >4%)
max_rb_pct   <- 4         # Max ribosomal % (excludes >4%)
max_doublet  <- 0.4       # Remove cells with doublet score >= 0.4
min_cells    <- 3         # Remove genes detected in < 3 cells
```

**LC-NE Marker Genes:**
- **Dbh** (Dopamine β-hydroxylase)
- **Th** (Tyrosine hydroxylase)
- **Slc6a2** (Norepinephrine transporter)
- **Slc18a2** (Vesicular monoamine transporter 2, VMAT2)

**Normalization:**
- Method: Counts-per-10k (CP10k) normalization followed by log-normalization
- Variable features: 2,000 genes (LC subset; note the Python/scVI transcriptomic step uses 1,500 HVGs)

**Batch Correction:**
- Algorithm: Harmony
- Batch variable: Sample identifier

## Environment

**R Version:** 4.4.2

**Key Dependencies:**
- Seurat 5.x (v5.4.0; SeuratObject v5.3.0)
- harmony
- tidyverse (dplyr, ggplot2)
- SeuratDisk or zellkonverter (for h5ad export)

**Resource Requirements:**
- Compute: Large instance (due to full dataset processing)
- Memory: ~16-32 GB recommended

## Output Structure

Results are written to `/results/`:

```
/results/
├── objects/
│   ├── LC_clustered_harmony.rds       # Full dataset, batch-corrected
│   ├── LC_subset.rds                   # Final LC-NE cells only
│   └── LC_subset.h5ad                  # AnnData export
└── figures/
    ├── umap_clusters.png
    ├── umap_markers.png
    ├── violin_markers.png
    └── qc_metrics.png
```

## Usage

### Reproducible Run
Execute the full pipeline via the `/code/run` script:
```bash
bash /code/run
```

The pipeline runs all 6 scripts sequentially with validation checks between steps.

### Interactive (Cloud Workstation)
Run scripts individually in R console or RStudio:
```r
source("code/01_build_full_object.R")
source("code/02_define_LC_population.R")
source("code/03_LC_subclustering.R")
source("code/02_b_publication_fig.R")
source("code/04_export_LC_to_h5ad.R")
source("code/05_retroseeq_convert.r")
```

**Note:** Scripts must be run in order, as each depends on outputs from previous steps.

## Notes

- **Memory Management:** The full dataset is processed in step 1; subsequent steps work on LC subset (~5-10% of cells)
- **Batch Effects:** Harmony correction assumes sample-level batch effects; adjust `group.by.vars` if needed
- **Marker Filtering:** LC-NE identification is stringent (clusters retained require >1 log CP10k expression for all markers); adjust thresholds if yield is too low/high
- **renv:** Package environment managed via `renv.lock` (auto-disabled for reproducible runs)

## Troubleshooting

**Low LC-NE Cell Yield:**
- Check marker gene expression thresholds in `02_define_LC_population.R`
- Verify marker genes are present in dataset (`Dbh`, `Th`, `Slc6a2`, `Slc18a2`)

**Memory Errors:**
- Increase compute instance size
- Reduce number of variable features used for PCA

**Harmony Convergence Issues:**
- Reduce `max.iter.harmony` parameter
- Check for batch variables with single samples

## Citation

If using this preprocessing pipeline, please cite:
- Seurat: Hao et al., Cell 2021
- Harmony: Korsunsky et al., Nature Methods 2019
- SeuratDisk/zellkonverter: For h5ad conversion

```
@article{xxxx,
  title   = {Topographic structure and function of locus coeruleus
norepinephrine neurons},
  author  = {Zhixiao Su},
  journal = {xxx},
  volume  = {xx},
  number  = {xx},
  pages   = {xxx},
  year    = {xxx},
  publisher = {xxx}
}
```
## Contact

For questions or issues, please open an issue on the [GitHub repository](https://github.com/AllenNeuralDynamics/LCNE_10x_preprocessing_R).