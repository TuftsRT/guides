---
tags: bioinformatics
---

# Set Up Conda Environment and Create Jupyter Kernel for scRNA-seq Analysis

Shirley Li, Bioinformatician, TTS Research Technology
xue.li37@tufts.edu

Date: 2024-11-01

## Overview

In this tutorial, you will learn how to:

- Create a Conda environment for single-cell RNA-seq analysis using Python-only packages.
- Install popular Python packages for scRNA-seq, such as Scanpy, and Scrublet.
- Set up a Jupyter kernel that uses the Conda environment for easy access to the tools in a notebook interface.

### Create a Conda Environment for scRNA-seq

1. Load `miniforge` and `conda-env-mod` module

```
module load miniforge/24.7.1-py312
module load conda-env-mod/default
```

2. Configure your conda

   **Note (steps in this session only needs to be executed ONCE)**

   Since your home directory has limited storage, it’s recommended to install conda packages in your group research storage space. Follow these steps:

   Create two directories in your group research storage space (one for storing the envs, one for storing the pkgs, for example: condaenv, condapkg)

   ```
   mkdir /cluster/tufts/XXXXlab/$USER/condaenv/
   mkdir /cluster/tufts/XXXXlab/$USER/condapkg/
   ```

   If you haven’t used conda before on the cluster, create a file named “.condarc” in your home directory.

   Now add the following 4 lines to the `.condarc` file in your home directory (modify according to your real path to the directories):

   ```
   envs_dirs:
     - /cluster/tufts/XXXXlab/$USER/condaenv/
   pkgs_dirs:
     - /cluster/tufts/XXXXlab/$USER/condapkg/
   ```

   After this, your `.condarc` file should look like this:

   ```
   envs_dirs:
     - /cluster/tufts/XXXXlab/$USER/condaenv/
   pkgs_dirs:
     - /cluster/tufts/XXXXlab/$USER/condapkg/
   channels:
     - bioconda
     - conda-forge
     - defaults
   ```

1. Create your conda environment with `conda-env-mod`

```
cd /cluster/tufts/XXXXlab/$USER/condaenv/
conda-env-mod create -p scrna_seq_py_env python=3.8  --jupyter
```

​ You will see something like this, and enter `y` to continue

```
  ...

The following NEW packages will be INSTALLED:

  _libgcc_mutex      conda-forge/linux-64::_libgcc_mutex-0.1-conda_forge
  _openmp_mutex      conda-forge/linux-64::_openmp_mutex-4.5-2_gnu
  asttokens          conda-forge/noarch::asttokens-2.4.1-pyhd8ed1ab_0
  bzip2              conda-forge/linux-64::bzip2-1.0.8-hd590300_5
  ca-certificates    conda-forge/linux-64::ca-certificates-2024.7.4-hbcca054_0
  ...

Proceed ([y]/n)? y
```

​ When it’s complete, you will see something like this:

```
...
Preparing transaction: ...working... done
Verifying transaction: ...working... done
Executing transaction: ...working... done
+---------------------------------------------------------------+
| To use this environment, load the following modules:          |
|     module load use.own                                       |
|     module load conda-env/scrna_seq_py_env-py3.12.5           |
| (then standard 'conda install' / 'pip install' / run scripts) |
+---------------------------------------------------------------+
```

### Install Selected Python Packages

1. Activate conda environment and install new packages

   ```
   module load use.own
   module load conda-env/scrna_seq_py_env-py3.12.5

   conda list # check packages installed in this environment

   pip install jupyter
   pip install numpy
   pip install pandas
   pip install anndata
   conda install -c conda-forge scanpy
   conda install -c bioconda scrublet
   pip install harmony-pytorch
   pip install gseapy
   pip install scanorama
   pip install pyscenic
   pip install scvi-tools
   pip install -i https://test.pypi.org/simple/ memento
   pip install pooch
   conda install -c conda-forge python-igraph

   conda list # check again
   ```

1. Create a jupyter kernel

   ```
   conda-env-mod kernel -n scrna_seq_py_env
   ```

   You will see something like this:

   ```
   requested kernel with arguments:  -n 'scrna_seq_py_env' --

   Jupyter kernel created: "Python (My scrna_seq_py_env Kernel)"
   +---------------------------------------------------------------+
   | We recommend installing packages into your kernel environment |
   | via the command line (with 'conda install' or 'pip install'). |
   +---------------------------------------------------------------+
   ```

## Using Open OnDemand Jupyter Lab

Natigate to [Open Ondemand](https://ondemand.pax.tufts.edu/)

In Open Ondemand dashboard, let’s go to `Interactive APPs` => `Jupyter Lab` and select the `number of hours`, `number of cores`, and `Amount of memory` that you would like to request and Launch this job.

Under `Notebook`, select the kernel you just created. Ex: `scrna_seq_py_env`

Start your python code from there.

Example code to check the installation:

```
# Import installed packages
import os
import seaborn as sns
import scanpy as sc
import scrublet as scr
import anndata
import harmony
import memento
import numpy as np
import pandas as pd
import scvi
import matplotlib.pyplot as plt
```

## Single-Cell RNA-seq Analysis Packages

### Scanpy

- **Summary**: `Scanpy` is a widely used Python package for analyzing large-scale single-cell RNA-seq datasets. It is optimized for scalability and supports workflows for preprocessing, clustering, dimensionality reduction, differential expression, and visualization of single-cell data.
- **Paper**: Wolf, F. A., Angerer, P., & Theis, F. J. (2018). "Scanpy: large-scale single-cell gene expression data analysis." _Genome Biology_, 19(1), 15. [https://doi.org/10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0) % codespell:ignore theis
- **Website**: [https://scanpy.readthedocs.io](https://scanpy.readthedocs.io)

### Scrublet

- **Summary**: `Scrublet` is a Python tool designed to detect doublets in single-cell RNA-seq data. Doublets are instances where two cells are captured in a single droplet, which can distort downstream analysis. Scrublet uses a k-nearest neighbors approach to identify and score potential doublets.
- **Paper**: Wolock, S. L., Lopez, R., & Klein, A. M. (2019). "Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data." _Cell Systems_, 8(4), 281–291.e9. [https://doi.org/10.1016/j.cels.2018.11.005](https://doi.org/10.1016/j.cels.2018.11.005)
- **Website**: [https://github.com/AllonKleinLab/scrublet](https://github.com/AllonKleinLab/scrublet)

### AnnData

- **Summary**: `AnnData` is a Python package that provides a framework for managing annotated data matrices, tailored for large-scale single-cell RNA-seq data. AnnData is widely used as the primary data structure in `Scanpy`, enabling efficient storage and handling of both raw and processed single-cell data.
- **Paper**: Virshup, I., Rybakov, S., Theis, F. J., Angerer, P., & Wolf, F. A. (2024). "anndata: Access and store annotated data matrices." _The Journal of Open Source Software_. [https://doi.org/10.21105/joss.04371](https://doi.org/10.21105/joss.04371) % codespell:ignore theis
- **Website**: [https://anndata.readthedocs.io](https://anndata.readthedocs.io)

### Harmony

- **Summary**: `Harmony` is a tool designed for batch effect correction in single-cell RNA-seq datasets. It integrates datasets from different batches or conditions by aligning data in a shared embedding space, allowing biological variation to be preserved while minimizing technical differences.
- **Paper**: Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F., Wei, K., ... & Raychaudhuri, S. (2019). "Fast, sensitive and accurate integration of single-cell data with Harmony." _Nature Methods_, 16(12), 1289-1296. [https://doi.org/10.1038/s41592-019-0619-0](https://doi.org/10.1038/s41592-019-0619-0)
- **Website**: [https://portals.broadinstitute.org/harmony](https://portals.broadinstitute.org/harmony)

### Memento

- **Summary**: `Memento` is a statistical tool tailored for single-cell RNA sequencing (scRNA-seq) analysis, with a focus on decoupling measurement noise from biological expression variability, thereby improving accuracy in differential expression studies.
- **Paper**: Kim, M. C., Gate, R., Lee, D. S., Marson, A., Ntranos, V., Ye, C. J. (2024). "Method of moments framework for differential expression analysis of single-cell RNA sequencing data." _Cell_, 187(22), P6393-6410.E16. [https://doi.org/10.1016/j.cell.2024.08.022](https://doi.org/10.1016/j.cell.2024.08.022)https://doi.org/10.1038/s41592-021-01125-y)
- **Website**: [https://github.com/yelabucsf/scrna-parameter-estimation](https://github.com/yelabucsf/scrna-parameter-estimation)

### scVI-tools

- **Summary**: `scVI-tools` is a framework built on top of PyTorch for scalable probabilistic modeling of single-cell data. It includes various models like scVI (single-cell variational inference), totalVI, and PEAKVI, used for data integration, dimensionality reduction, differential expression, and multi-omics data analysis.
- **Paper**: Gayoso, A., Lopez, R., Xing, G., Boyeau, P., Wu, K., Jayasuriya, M., et al. (2022). "A Python library for probabilistic analysis of single-cell omics data." _Nature Biotechnology_, 40, 163–166. [https://www.nature.com/articles/s41587-021-01206-w](https://www.nature.com/articles/s41587-021-01206-w)
- **Website**: [https://scvi-tools.org](https://scvi-tools.org)
