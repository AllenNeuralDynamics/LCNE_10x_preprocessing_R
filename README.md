# LCNE 10x snRNA-seq Preprocessing Pipeline

**Locus Coeruleus Norepinephrine (LC-NE) Cell Identification and Analysis**

This capsule performs sequential preprocessing of 10x Genomics single-nucleus RNA-seq data to identify, filter, and characterize locus coeruleus norepinephrine (LC-NE) neurons using Seurat 5 and Harmony batch correction.

## Overview

The pipeline implements a multi-stage workflow:
1. Quality control and Harmony batch correction on full dataset
2. LC-NE population identification via marker gene expression
3. LC-specific subclustering and outlier removal
4. Publication-quality visualization
5. Export to AnnData (h5ad) format for Python workflows

## Input Data

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
- Quality control filtering (gene counts, mitochondrial %, ribosomal %)
- SCTransform normalization
- Harmony batch correction across samples
- Initial dimensionality reduction (PCA, UMAP)
- Preliminary clustering

**Output:** `LC_clustered_harmony.rds` (Seurat object) -> this creates ~ 6GB of data. 

### Script 2: Define LC Population (`02_define_LC_population.R`)
**Purpose:** Identify LC-NE neurons based on marker genes

**Steps:**
- Load batch-corrected Seurat object
- Screen clusters for LC-NE markers (Dbh, Th, Slc6a2)
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

## Key Parameters

**Quality Control:**
```r
min_genes <- 200          # Minimum genes per cell
max_genes <- 5000         # Maximum genes per cell
max_mt_pct <- 10          # Max mitochondrial % (excludes >10%)
max_rb_pct <- 50          # Max ribosomal %
```

**LC-NE Marker Genes:**
- **Dbh** (Dopamine β-hydroxylase)
- **Th** (Tyrosine hydroxylase)
- **Slc6a2** (Norepinephrine transporter)

**Normalization:**
- Method: SCTransform (Seurat v5)
- Variable features: 3,000 genes

**Batch Correction:**
- Algorithm: Harmony
- Batch variable: Sample identifier

## Environment

**R Version:** 4.4.2

**Key Dependencies:**
- Seurat 5.x
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

The pipeline runs all 5 scripts sequentially with validation checks between steps.

### Interactive (Cloud Workstation)
Run scripts individually in R console or RStudio:
```r
source("code/01_build_full_object.R")
source("code/02_define_LC_population.R")
source("code/03_LC_subclustering.R")
source("code/02_b_publication_fig.R")
source("code/04_export_LC_to_h5ad.R")
```

**Note:** Scripts must be run in order, as each depends on outputs from previous steps.

## Notes

- **Memory Management:** The full dataset is processed in step 1; subsequent steps work on LC subset (~5-10% of cells)
- **Batch Effects:** Harmony correction assumes sample-level batch effects; adjust `group.by.vars` if needed
- **Marker Filtering:** LC-NE identification is stringent; adjust thresholds if yield is too low/high
- **renv:** Package environment managed via `renv.lock` (auto-disabled for reproducible runs)

## Troubleshooting

**Low LC-NE Cell Yield:**
- Check marker gene expression thresholds in `02_define_LC_population.R`
- Verify marker genes are present in dataset (`Dbh`, `Th`, `Slc6a2`)

**Memory Errors:**
- Increase compute instance size
- Reduce number of variable features in SCTransform

**Harmony Convergence Issues:**
- Reduce `max.iter.harmony` parameter
- Check for batch variables with single samples

## Citation

If using this preprocessing pipeline, please cite:
- Seurat: Hao et al., Cell 2021
- Harmony: Korsunsky et al., Nature Methods 2019
- SeuratDisk/zellkonverter: For h5ad conversion

## Contact

For questions or issues, please open an issue on the [GitHub repository](https://github.com/AllenNeuralDynamics/LCNE_10x_preprocessing_R).