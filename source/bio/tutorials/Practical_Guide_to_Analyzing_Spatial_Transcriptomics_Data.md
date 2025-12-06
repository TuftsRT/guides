# Practical Guide to Analyzing Spatial Transcriptomics Data

Author: Shirley Li, xue.li37@tufts.edu               
Date: 2025-12-06

## 1. Overview

This document introduces three major spatial transcriptomics platforms used in modern biological and biomedical research:

- **Vizgen MERFISH**: Single-molecule, imaging-based, single-cell resolved
- **10x Genomics Visium**: Sequencing-based, spot-level whole-transcriptome
- **NanoString GeoMx DSP**: Region-of-interest (ROI) based spatial profiling

The focus is on **practical application**, including:

- What each platform is used for
- What types of files they generate
- How to analyze the data
- Fully reproducible **R-based workflows**

This is not a theory-focused document; instead, it is designed to help researchers analyze data on the Tufts HPC environment or locally.


# 2. Vizgen MERFISH

## 2.1 What MERFISH Is

MERFISH (Multiplexed Error-Robust FISH) from Vizgen is an **imaging-based single-cell spatial transcriptomics** technology. It detects individual RNA molecules directly in tissue, enabling subcellular precision.

Key characteristics:

- Single-cell resolution
- High spatial fidelity
- Molecule-level detection
- Customizable gene panels
- Produces cell-by-gene matrices and molecule coordinates



## 2.2 Typical Output Files

```
cell_by_gene.csv
cell_metadata.csv
molecule_list.csv
fov_positions.json
images/
```

The primary data for analysis:

- `cell_by_gene.csv`: counts matrix
- `cell_metadata.csv`: x/y positions, cell segmentation, QC metrics



## 2.3 MERFISH Analysis Workflow in R

### Load required libraries

```
library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)
```

### Step 1. Load counts and metadata

```
counts <- read_csv("cell_by_gene.csv")
counts <- as.data.frame(counts)
rownames(counts) <- counts$cell
counts$cell <- NULL
counts <- t(counts)

metadata <- read_csv("cell_metadata.csv")
rownames(metadata) <- metadata$cell

obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
```

### Step 2. QC and filtering

```
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

obj <- subset(
    obj,
    subset = nFeature_RNA > 30 &
             nCount_RNA > 500 &
             percent.mt < 25
)
```

### Step 3. Normalization and scaling

```
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
```

### Step 4. PCA, neighbors, clustering

```
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj)
obj <- RunUMAP(obj, dims = 1:20)
DimPlot(obj, reduction = "umap")
```

### Step 5. Spatial visualization

Requires `x` and `y` coordinates in metadata.

```
ggplot(obj@meta.data, aes(x = x, y = y, color = seurat_clusters)) +
    geom_point(size = 0.6) +
    coord_fixed() +
    theme_bw()
```

### Step 6. Marker genes and annotation

```
markers <- FindAllMarkers(obj)
```

### Step 7. Neighborhood analysis

Example using nearest neighbors.

```
library(FNN)
coords <- obj@meta.data[, c("x", "y")]
nn <- get.knn(coords, k = 10)
```



# 3. 10x Genomics Visium

## 3.1 What Visium Is

10x Visium is a **sequencing-based spatial transcriptomics** system. mRNA is captured on spatially barcoded spots (~55 µm), processed through Space Ranger, and analyzed via R or Python.

Key characteristics:

- Whole-transcriptome profiling
- Multiple cells per spot
- Supports FFPE and fresh frozen samples
- Strong for tissue-level discovery and regional differences



## 3.2 Output Files from Space Ranger

```
filtered_feature_bc_matrix/
spatial/
analysis/
```

Files used for analysis:

- `matrix.mtx.gz`
- `features.tsv.gz`
- `barcodes.tsv.gz`
- `tissue_positions_list.csv`
- high-resolution H&E image



## 3.3 Visium Analysis Workflow in R

### Step 1. Load data

```
library(Seurat)

obj <- Load10X_Spatial(data.dir = "path/to/visium_folder")
```

### Step 2. QC

```
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

obj <- subset(
    obj,
    subset = nCount_Spatial > 500 &
             nFeature_Spatial > 200 &
             percent.mt < 25
)
```

### Step 3. Normalization and clustering

```
obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj)
obj <- RunUMAP(obj, dims = 1:20)

SpatialDimPlot(obj, label = TRUE)
```

### Step 4. Differential expression

```
markers <- FindAllMarkers(obj)
```

### Step 5. Spatially variable features

```
svg <- FindSpatiallyVariableFeatures(obj, assay = "SCT")
```

### Step 6. Cell type deconvolution (example: SPOTlight)

```
library(SPOTlight)

# sc is a Seurat object with scRNA-seq reference
res <- spotlight_deconvolution(se_sc = sc, st_obj = obj)
```

### Step 7. Pathway enrichment

```
library(clusterProfiler)
library(org.Hs.eg.db)

genes <- markers$gene[markers$cluster == "0"]

enrich <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL"
)
```



# 4. NanoString GeoMx DSP

## 4.1 What GeoMx DSP Is

GeoMx DSP is an ROI-based spatial profiling technology used widely in translational and pathology settings. ROIs are defined using immunofluorescence, and gene expression is collected via UV cleavage of barcoded probes.

Key characteristics:

- ROI-level profiling, not single-cell
- Works well with FFPE
- Supports whole transcriptome and targeted panels
- Ideal for comparisons between regions, tumor microenvironment, inflamed vs non-inflamed areas



## 4.2 Typical GeoMx Files

```
exprMat.txt
metadata.txt
panel_annotation.txt
QC_stats/
image_tiles/
```



## 4.3 GeoMx Analysis Workflow in R

### Step 1. Load data

```
library(GeoMxTools)

obj <- readNanoStringGeoMxSet(
    dccFiles = "path/to/dcc/",
    pkcFiles = "path/to/pkc/",
    phenoDataFile = "metadata.txt",
    exprsFile = "exprMat.txt"
)
```

### Step 2. QC

```
obj <- shiftCountsOne(obj)
obj <- setSegmentQCFlags(obj)
obj <- setGeneQCFlags(obj)

obj <- obj[, obj$SegmentQC == "PASS"]
```

### Step 3. Normalization

```
obj <- normalize(obj, norm_method = "quant")
```

### Step 4. PCA

```
obj <- runPCA(obj)
plotPCA(obj)
```

### Step 5. Differential expression via limma

```
library(limma)

design <- model.matrix(~ group, data = pData(obj))
fit <- lmFit(exprs(obj), design)
fit <- eBayes(fit)

topTable(fit, coef = 2)
```

### Step 6. Pathway analysis

```
library(clusterProfiler)
genes <- rownames(topTable(fit, coef = 2, number = Inf, p.value = 0.05))

enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL"
)
```



# 5. Summary of Recommended Tools

| Platform  | Resolution              | Primary Tools (R)                  | Typical Analyses                                          |
| --------- | ----------------------- | ---------------------------------- | --------------------------------------------------------- |
| MERFISH   | Single-cell/subcellular | Seurat, ggplot2, FNN               | Cell classification, spatial visualization, neighborhoods |
| Visium    | Spot-level              | Seurat, SPOTlight, clusterProfiler | Clustering, deconvolution, spatially variable genes       |
| GeoMx DSP | ROI-level               | GeoMxTools, limma, clusterProfiler | Region comparisons, FFPE profiling, pathway analysis      |



# 6. Notes                  
- All provided code is compatible with SLURM workloads on the Tufts HPC.
- Use the latest Open OnDemand RStudio session for interactive development.
- We recommend using Seurat v5 for spatial transcriptomics pipelines.
- For MERFISH or Visium analyses, request 16–32 CPUs as a start. Always pair with sufficient memory (--mem=64G or more).
- Data import and QC can be time-consuming. Save your Seurat (or other) object to an .rds file and continue from the saved object to avoid re-running preprocessing.
