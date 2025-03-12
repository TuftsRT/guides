# A Beginner’s Guide to Differential Expression Analysis with DESeq2 in R on Open OnDemand

## Prerequisite

1. **Access to Tufts HPC**: If you don’t already have an account, you can request one [here](https://shib-idp.tufts.edu/idp/profile/SAML2/POST/SSO?execution=e1s1).
1. **Lab Storage Access**: For analysis, we recommend using your lab’s designated storage rather than the `$HOME` directory. For more details on storage options, refer to [this guide](https://it.tufts.edu/research-data-storage).

## Reference

- [**DESeq2 official website**](https://bioconductor.org/packages/release/bioc/html/DESeq2.html): Comprehensive information on DESeq2, including installation instructions, documentation, and updates.

- [User guide for how to use Tufts HPC](https://tufts.app.box.com/v/Pax-User-Guide): Detailed instructions for accessing and using Tufts High-Performance Computing (HPC) resources.

- DESeq2 tutorials

  - [A beginner-friendly guide to using DESeq2 for differential gene expression analysis. ](https://lashlock.github.io/compbio/R_presentation.html)
  - [Step-by-step walkthrough for DESeq2 analysis. ](https://introtogenomics.readthedocs.io/en/latest/2021.11.11.DeseqTutorial.html)
  - [Differential expression with DEseq2. ](https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/)
  - [Harvard Chan Bioinformatics Core Training: Introduction to DGE. ](https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html)

- [How to Download and Install R and RStudio](https://rstudio-education.github.io/hopr/starting.html)

  - If you are using your own laptop, follow the instructions in this tutorial to install R and RStudio.

- [Install Bioconductor Packages](https://www.bioconductor.org/install/)

  - If you are not using our Open OnDemand RStudio_Pax app, follow the instructions in this tutorial to install `DESeq2` packages.

## Step1: Access the Cluster via Open OnDemand

**Log in** to [Open OnDemand](https://ondemand.pax.tufts.edu) using your Tufts HPC credentials.

**Start an RStudio Job**:

1. Navigate to **Interactive Apps**, then select **RStudio_Pax**.

1. Fill in the required fields:

   - **Number of hours**: Specify the time needed for your session. (Recommended starting point: 24 hours)
   - **Number of CPU cores**: Choose based on your computational needs. (Recommended starting point: 8 cores)
   - **Amount of memory**: Allocate sufficient memory for your analysis. (Recommended starting point: 32 GB)
   - **R Version**: Select `4.4.1`.

1. Click **Launch** to start the job.

**Note: guidelines for selecting computational resources (cores and memory) based on dataset size when performing DESeq2 differential expression analysis.**

- **Cores**:

  - Small datasets (\<10 samples): 1–2 cores.

  - Medium datasets (10–50 samples): 4–8 cores.

  - Large datasets (>50 samples): 8–16 cores, especially if using parallelization (e.g., `BiocParallel`).

- **Memory**:

  - Small datasets: 8 GB.

  - Medium datasets: 16–32 GB.

  - Large datasets: 64 GB or more for large matrices or numerous samples.

## Step2: Load the DESeq2 package and other necessary libraries

The version of DESeq2 currently pre-installed in R is **1.44.0**.

Verify it using the following commands in R:

```
library(DESeq2)
packageVersion("DESeq2")
```

If you need a newer version of DESeq2 (e.g., 1.46.0), you may need to update the package manually or request installation by the system administrator.

Load other libraries:

```
library(DESeq2)
library(ggplot2)
library(pheatmap)
```

## Step3: Begin the Analysis

### 3.1 Input Data Preparation

#### Requirements:

- A **count matrix** (rows = genes, columns = samples).
- A **metadata file** with sample information (e.g., condition, treatment).

#### Example:

Count Matrix (`counts_matrix.csv`):

| GeneID | Sample1 | Sample2 | Sample3 | Sample4 |
| ------ | ------- | ------- | ------- | ------- |
| GeneA  | 100     | 120     | 90      | 80      |
| GeneB  | 200     | 190     | 210     | 205     |

Metadata (`metadata.csv`):

| Sample  | Condition |
| ------- | --------- |
| Sample1 | Control   |
| Sample2 | Control   |
| Sample3 | Treatment |
| Sample4 | Treatment |

#### Load Data:

```
# Load count matrix
counts <- read.csv("counts_matrix.csv", row.names = 1)

# Load metadata
metadata <- read.csv("metadata.csv", row.names = 1)

# Ensure row names of metadata match column names of counts
metadata <- metadata[match(colnames(counts), rownames(metadata)), ]
```

### 3.2 Create DESeq2 Dataset

Create the DESeq2 object using `DESeqDataSetFromMatrix`:

```
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ Condition
)

# Pre-filtering: Remove low count genes
dds <- dds[rowSums(counts(dds)) > 10, ]
```

### 3.3 Differential Expression Analysis

```
dds <- DESeq(dds)
results <- results(dds)
results <- results[order(results$padj), ]
head(results)
```

### 3.4 Exploratory Data Analysis

```
# PCA plot
vsd <- vst(dds, blind = TRUE)  # Variance-stabilizing transformation
plotPCA(vsd, intgroup = "Condition")

# Heatmap of Sample Distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix, main = "Sample Distance Heatmap")
```

### 3.5 Volcano Plot

```
# Add significance column
results$significant <- results$padj < 0.05 & abs(results$log2FoldChange) > 1

# Volcano plot
ggplot(results, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10 P-value")
```

### 3.6 Annotating Results

```
# Add gene symbols or annotations to results (if available):
annotated_results <- merge(results, gene_annotations, by = "GeneID", all.x = TRUE)

# Save results to a CSV file:
write.csv(annotated_results, "DESeq2_results.csv")
```
