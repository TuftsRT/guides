---
tags: bioinformatics
---

# scRNA-seq analysis with RStudio

Shirley Li, Bioinformatician, TTS Research Technology
xue.li37@tufts.edu

Date: 2024-11-01

## Overview

In this tutorial, you will learn how to:

- Set up an R environment for single-cell RNA-seq analysis.
- Install and configure popular R packages for scRNA-seq, such as **Seurat**, **SingleCellExperiment**, and **scater**.
- Load, preprocess, and analyze scRNA-seq data in R with Seurat.

## Set Up RStudio on Open OnDemand for scRNA-seq

1. **Open OnDemand RStudio App**

   - Log in to [Open OnDemand](https://ondemand.pax.tufts.edu/) with your UTLN
   - You will see `RStudio Pax` under `Interactive Apps` and other topic-specific Rstudio apps under `Bioinformatics Apps`
   - Fill in the parameters according to your needs. Start with 64GB and 12 cores.
   - Launch the job.

1. **Load Packages in R**:

   ```r
   # Load the necessary libraries for scRNA-seq analysis
   library(Seurat)
   library(SeuratData)
   library(SingleCellExperiment)
   library(scater)
   library(scran)
   # code to install scRNAseq
   # BiocManager::install("scRNAseq")
   library(scRNAseq)
   library(monocle)
   ```

## Single-Cell RNA-seq Analysis Packages

### Seurat

- **Summary**: `Seurat` is a widely used R package for single-cell RNA-seq analysis. It supports preprocessing, clustering, dimensionality reduction, differential expression, and visualization of single-cell data.
- **Paper**: Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck III, W. M., et al. (2019). "Comprehensive Integration of Single-Cell Data." *Cell*, 177(7), 1888–1902. [https://doi.org/10.1016/j.cell.2019.05.031](https://doi.org/10.1016/j.cell.2019.05.031)
- **Website**: [https://satijalab.org/seurat/](https://satijalab.org/seurat/)

### SingleCellExperiment

- **Summary**: `SingleCellExperiment` provides a flexible framework for representing single-cell data in R, including assays, row/column metadata, and reduced dimensions. It serves as a foundational data structure in the Bioconductor ecosystem.
- **Paper**: Lun, A. T. L., & Risso, D. (2023). "SingleCellExperiment: an object class for single-cell analysis in R." *F1000Research*, 7, 521. [https://doi.org/10.12688/f1000research.16298.2](https://doi.org/10.12688/f1000research.16298.2)
- **Website**: [https://bioconductor.org/packages/SingleCellExperiment/](https://bioconductor.org/packages/SingleCellExperiment/)

### scater

- **Summary**: `scater` is used for pre-processing and quality control of single-cell data. It helps to generate various visualizations and metrics for understanding the quality and variability within scRNA-seq datasets.
- **Paper**: McCarthy, D. J., Campbell, K. R., Lun, A. T. L., & Wills, Q. F. (2017). "Scater: Pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R." *Bioinformatics*, 33(8), 1179–1186. [https://doi.org/10.1093/bioinformatics/btw777](https://doi.org/10.1093/bioinformatics/btw777)
- **Website**: [https://bioconductor.org/packages/scater/](https://bioconductor.org/packages/scater/)

### scran

- **Summary**: `scran` provides efficient methods for single-cell RNA-seq data normalization, clustering, and marker detection. It’s optimized for scalability with large datasets and includes several statistical tools for single-cell analysis.
- **Paper**: Lun, A. T. L., & Marioni, J. C. (2017). "Overcoming systematic errors caused by log-transformation of normalized single-cell RNA sequencing data." *Nature Methods*, 13, 795–798. [https://doi.org/10.1038/nmeth.3963](https://doi.org/10.1038/nmeth.3963)
- **Website**: [https://bioconductor.org/packages/scran/](https://bioconductor.org/packages/scran/)

### scRNAseq

- **Summary**: `scRNAseq` is a data package that includes several example scRNA-seq datasets for practice and benchmarking of analytical pipelines.
- **Paper**: Available on Bioconductor.
- **Website**: [https://bioconductor.org/packages/scRNAseq/](https://bioconductor.org/packages/scRNAseq/)

### Monocle

- **Summary**: `Monocle` is a package designed for analyzing single-cell trajectories. It identifies and orders cells based on gene expression dynamics over pseudotime, which can help reveal cell differentiation and developmental processes.

- **Paper**: Trapnell, C., Cacchiarelli, D., Grimsby, J., Pokharel, P., Li, S., Morse, M., et al. (2014). "The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells." *Nature Biotechnology*, 32(4), 381–386. https://doi.org/10.1038/nbt.2859

- **Website**: http://cole-trapnell-lab.github.io/monocle-release/

## Example Code for Data Loading and Visualization

Below is a step-by-step example code for loading, preprocessing, and visualizing single-cell RNA-seq data with Seurat, including comments to explain each step:

```r
# Load packages
library(Seurat)
library(SeuratData)

# Load the pbmc3k dataset (replace with your dataset if needed)
SeuratData::InstallData("pbmc3k")
data("pbmc3k") # Replace with your dataset

# Create a Seurat object
pbmc_counts <- as.matrix(pbmc3k@assays$RNA@counts)
pbmc <- CreateSeuratObject(counts = pbmc_counts)

# Quality control and filtering
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

# Normalize the data
# This step adjusts for differences in sequencing depth across cells.
pbmc <- NormalizeData(pbmc)

# Identify variable features
# Variable features are genes that exhibit high variability across cells.
# These features are used in downstream analyses to focus on informative genes.
pbmc <- FindVariableFeatures(pbmc)

# Scale the data
# Scaling centers and scales each gene, making them comparable for PCA and clustering.
# This is an important step for most dimensionality reduction techniques.
pbmc <- ScaleData(pbmc)

# Step 4: Run Principal Component Analysis (PCA)
# PCA reduces the dimensionality of the data, allowing us to identify major sources of variation.
# Here, we use the variable features identified earlier to perform PCA.
pbmc <- RunPCA(pbmc)

# Step 5: Find Neighbors
# This step identifies nearest neighbors for each cell based on their PCA scores.
# It is an essential step before clustering the cells.
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Step 6: Cluster the cells
# Clustering groups cells with similar expression profiles, aiding in cell-type identification.
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Step 7: Run UMAP for visualization
# UMAP (Uniform Manifold Approximation and Projection) is a popular method for visualizing high-dimensional data.
# This step reduces the data to two dimensions, making it easier to visualize clusters.
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Step 8: Plot the UMAP results
# This visualization shows the clusters identified by Seurat, each in a different color.
DimPlot(pbmc, reduction = "umap")
```

This annotated code provides a guide for basic steps of single-cell RNA-seq analysis in Seurat, from data normalization to visualization.

**For a detailed tutorial on using Seurat, please visit [this link](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)**

## Memory and Core Requirements

### **Memory Requirements**

- **Small Dataset (up to ~10,000 cells)**: 16–32 GB RAM should suffice for typical preprocessing and clustering tasks.
- **Medium Dataset (up to ~50,000 cells)**: 32–64 GB RAM is recommended to handle most analyses comfortably.
- **Large Dataset (over ~100,000 cells)**: 128 GB or more might be needed, particularly for steps like dimensionality reduction and integration.

### **Core Requirements**

- **Small Analyses**: For small datasets or exploratory analysis, 2–4 cores are generally sufficient.
- **Medium to Large Analyses**: For larger datasets, 8–16 cores are beneficial, especially for parallelizable steps such as normalization, dimensionality reduction, and clustering.
- **Very Large or High-throughput Analyses**: Using 24 cores or more can help expedite analyses on large datasets, especially when utilizing parallelized workflows with packages like `BiocParallel`.

### Example Setup Recommendations

For a typical medium-sized scRNA-seq analysis in R, a setup with **64 GB of memory and 8–12 cores** is often effective. If you’re working with very large datasets or performing computationally intensive steps like integration across multiple datasets, consider **128 GB of memory and 16 or more cores** to optimize performance and speed.
