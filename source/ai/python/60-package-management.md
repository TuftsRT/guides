---
tags: python conda hpc
---

# Python Package Management with Conda and Pip

Managing Python packages and environments is one of the most important skills for reproducible research. This guide covers the tools and workflows you need to keep your Python projects organized, stable, and shareable.

## Why Package Management Matters

Python packages are constantly updated, and different projects may require different versions of the same package. Without careful environment management:

- Updating a package for one project can break another
- Code that works on your machine may not work on a collaborator's machine or on the HPC
- Reproducing old results months or years later becomes difficult

The solution is to create **isolated environments**, one per project, each containing a specific, recorded set of package versions.

---

## Conda vs. Pip: Which Should I Use?

Both conda and pip install Python packages, but they work differently:

|                            | conda                                                         | pip                                                  |
| -------------------------- | ------------------------------------------------------------- | ---------------------------------------------------- |
| **What it installs**       | Python packages + non-Python dependencies (C libraries, etc.) | Python packages only                                 |
| **Environment management** | Built-in (`conda create`, `conda activate`)                   | Requires separate tool (venv or virtualenv)          |
| **Package repository**     | conda-forge, defaults                                         | PyPI                                                 |
| **Best for**               | Scientific/data science packages (NumPy, PyTorch, etc.)       | Packages not on conda channels; pure Python packages |
| **Recommended for HPC**    | Yes                                                           | Supplementary use within a conda env                 |

**Our recommendation:** Use conda as your primary package manager. When a package is not available on conda channels, you can still install it with pip *inside* an active conda environment.

> See [Which Python Setup is Right for You?](10-which-python-setup.md) for more context on why we recommend conda.

---

## Managing Conda Environments

### Create a New Environment

Create a named environment with a specific Python version:

```bash
conda create --name myproject python=3.11
```

Replace `myproject` with a descriptive name for your project. Environment names should not contain spaces.

### Activate and Deactivate

```bash
conda activate myproject    # Switch to this environment
conda deactivate            # Return to the base environment
```

### List All Environments

```bash
conda env list
```

### Install Packages

Always activate the environment first, then install:

```bash
conda activate myproject
conda install numpy pandas matplotlib scikit-learn
```

Install a specific version:

```bash
conda install numpy=1.26.4
```

Install from the conda-forge channel (wider selection of packages):

```bash
conda install -c conda-forge xarray
```

### Install with pip (when needed)

If a package is not available via conda, install it with pip while the conda environment is active:

```bash
conda activate myproject
pip install some-package
```

> **Important:** Always install as much as possible with conda before resorting to pip. Mixing conda and pip can sometimes cause conflicts. When you do use pip inside a conda environment, use it last, after all conda installs are done.

### Update and Remove Packages

```bash
conda update numpy           # Update a specific package
conda remove pandas          # Remove a package
```

### Delete an Environment

```bash
conda env remove --name myproject
```

---

## Sharing and Reproducing Environments

The key to reproducibility is recording exactly which packages and versions your project uses, so that someone else (or your future self) can recreate the same environment.

### Export an Environment

Export your current environment to a YAML file:

```bash
conda env export > environment.yaml
```

This creates a file listing all packages and their exact versions.

### Recreate an Environment from a YAML File

```bash
conda env create --file environment.yaml
```

> **Tip:** Add `environment.yaml` to your project's Git repository so collaborators can reproduce your environment with one command.

### Cross-Platform Environments

A full `conda env export` includes platform-specific build strings that may not work on different operating systems. For a more portable export, use:

```bash
conda env export --no-builds > environment.yaml
```

Or for a minimal specification that only lists the packages you explicitly installed (leaving conda to resolve compatible versions on any platform):

```bash
conda env export --from-history > environment.yaml
```

---

## Managing pip Requirements Files

If you are working in a context where pip is the primary tool (e.g., a project with a `requirements.txt` file), here is how to work with it effectively.

### Create a requirements.txt

In an active environment, export all installed packages:

```bash
pip freeze > requirements.txt
```

This creates a file like:

```
numpy==1.26.4
pandas==2.2.1
scikit-learn==1.4.0
```

### Install from requirements.txt

```bash
pip install -r requirements.txt
```

> **Tip:** When possible, use conda environments (with `environment.yaml`) instead of `requirements.txt` for data science projects, since conda handles non-Python dependencies more reliably.

---

## Conda on the Tufts HPC

Conda environments you create on the HPC are stored in your home directory and are available across all HPC sessions. To get started:

1. Connect to the HPC (via SSH or [OnDemand](https://ondemand-p01.pax.tufts.edu/))
1. Load the conda module if it is not already available:
   ```bash
   module load miniforge/25.3.0
   ```
   (Check available versions with `module av miniforge`)
1. Create and manage environments as described above

> **Home directory quota:** Conda environments can be large. If you have access to lab research storage (`/cluster/tufts/XXXXlab/`), configure conda to store environments and packages there to avoid filling your home directory quota. See the [HPC conda documentation](../../hpc/application/10-condaenv.md) for instructions.

To use your conda environment in a SLURM batch job, activate it in your job script:

```bash
#!/bin/bash
#SBATCH -p batch
#SBATCH -n 4
#SBATCH --mem=16g
#SBATCH -t 4:00:00

module purge
module load miniforge/25.3.0
source activate myproject

python my_analysis.py
```

For more information on conda on the HPC, see the [HPC conda documentation](../../hpc/application/10-condaenv.md).

---

## Getting Help

For package management questions, contact **datalab-support@elist.tufts.edu**.
