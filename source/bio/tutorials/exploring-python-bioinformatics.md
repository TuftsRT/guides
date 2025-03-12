---
tags: bioinformatics
---

# Exploring Python Bioinformatics Packages with Jupyter Notebook

Shirley Li, Bioinformatician, TTS Research Technology
xue.li37@tufts.edu

Date: 2024-11-01

In this tutorial, we will use the **Anndata** package as an example to show how to run interactive Python sessions through the Tufts Open OnDemand Jupyter Notebook.

## Prerequisite

1. Familiarity with Linux commands
1. Experience working with conda environments

## Creating conda environment

1. Start an interactive job session\
   `srun -p interactive -n 1 --time=04:00:00 --mem 4g --pty bash`

1. Load anacoda or minoconda module

   `module load anaconda/2021.05`

   or

   `module load anaconda/2021.11`

   or

   `module load miniconda/23.10`

1. Load conda-env-mod module

   `module load conda-env-mod/default`

1. Configure your conda

   **_NOTE (steps in this session only needs to be executed ONCE)_**

   Since your home directory has limited storage, it’s recommended to install conda packages in your group research storage space. Follow these steps:

   Create two directories in your group research storage space (one for storing the envs, one for storing the pkgs, for example: condaenv, condapkg)

   `$ mkdir /cluster/tufts/XXXXlab/$USER/condaenv/`

   `$ mkdir /cluster/tufts/XXXXlab/$USER/condapkg/`

   If you haven't used conda before on the cluster, create a file named ".condarc" in your home directory.

   Now add the following 4 lines to the `.condarc` file in your home directory (modify according to your real path to the directories):

   ```
   envs_dirs:
     - /cluster/tufts/XXXXlab/$USER/condaenv/
   pkgs_dirs:
     - /cluster/tufts/XXXXlab/$USER/condapkg/
   ```

   After this, your `.condarc` file should look like this

   `$ cat ~/.condarc`

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

1. Create your conda environment with conda-env-mod

   Change `yourenvname` to the name of the environment you intend to create

   ```
   cd /cluster/tufts/XXXXlab/$USER/condaenv/
   conda-env-mod create -p yourenvname python=3.8  --jupyter
   ```

   You will see something like this, and enter `y` to continue

   ```
   The following NEW packages will be INSTALLED:

     _libgcc_mutex      conda-forge/linux-64::_libgcc_mutex-0.1-conda_forge
     _openmp_mutex      conda-forge/linux-64::_openmp_mutex-4.5-2_gnu
     asttokens          conda-forge/noarch::asttokens-2.4.1-pyhd8ed1ab_0
     bzip2              conda-forge/linux-64::bzip2-1.0.8-hd590300_5
     ca-certificates    conda-forge/linux-64::ca-certificates-2024.7.4-hbcca054_0
     ...

   Proceed ([y]/n)? y
   ```

   When it's complete, you will see something like this.

   ```


   Preparing transaction: ...working... done
   Verifying transaction: ...working... done
   Executing transaction: ...working... done
   +---------------------------------------------------------------+
   | To use this environment, load the following modules:          |
   |     module load use.own                                       |
   |     module load conda-env/bio_test-py3.11.5                   |
   | (then standard 'conda install' / 'pip install' / run scripts) |
   +---------------------------------------------------------------+
   ```

1. Activate conda environment and install new packages

   Note: `conda-env/bio_test-py3.11.5 ` this may be different and it depends on what `yourenvname` you have

   ```
   module load use.own
   module load conda-env/bio_test-py3.11.5

   conda list # check packages installed in this environment

   pip install jupyter
   pip install anndata

   conda list # check again
   ```

​

```
# packages in environment at bio_test:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                       2_gnu    conda-forge
anndata                   0.10.8                   pypi_0    pypi
array-api-compat          1.7.1                    pypi_0    pypi
asttokens                 2.4.1              pyhd8ed1ab_0    conda-forge
bzip2                     1.0.8                hd590300_5    conda-forge
...
```

7. Create a jupyter kernel

   `conda-env-mod kernel -n bio_test`

You will see something like this:

```
Setting CONDA_ENVS_PATH=/cluster/home/xli37/.conda/envs/rhel7.8/conda-23.10.0
New environments will be created in this directory unless --prefix is specified.
requested kernel with arguments:  -n 'bio_test' --

Jupyter kernel created: "Python (My bio_test Kernel)"
+---------------------------------------------------------------+
| We recommend installing packages into your kernel environment |
| via the command line (with 'conda install' or 'pip install'). |
```

​

## Using Open Ondemand Jupyter Lab

Natigate to [Open Ondemand](https://ondemand.pax.tufts.edu/)

In Open Ondemand dashboard, let's go to `Interactive APPs` => `Jupyter Lab` and select the `number of hours`, `number of cores`, and `Amount of memory` that you would like to request and Launch this job.

Under `Notebook`, select the kernel you just created. Ex: anndata_python.

Start your python code from there.

Example code to check the Anndata installation:

```
import anndata as ad
from scipy.sparse import csr_matrix
print(ad.__version__)
```

## Tutorials for ANNDATA

https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html

## Some basic python commands

Check current path

```
import os
print(os.getcwd())
```

Go to a new path

```
os.chdir('/cluster/tufts/XXLAB/$USER/')
```

Check what files exist in current path

```
os.listdir(os.getcwd())
```
