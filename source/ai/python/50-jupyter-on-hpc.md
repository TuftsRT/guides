---
tags: python jupyter hpc conda
---

# Jupyter on the Tufts HPC Cluster

The Tufts High-Performance Computing (HPC) cluster supports Jupyter Lab through the **OnDemand** web portal. This lets you run Jupyter notebooks with access to the cluster's CPUs, memory, and GPUs, without needing to use the command line for setup.

This guide covers:

1. [Launching a Jupyter Lab session on the HPC](#launching-jupyter-lab-on-the-hpc)
1. [Selecting and managing Python kernels](#selecting-a-kernel)
1. [Using conda environments as Jupyter kernels](#using-a-conda-environment-as-a-kernel)
1. [Transferring files between your computer and the HPC](#transferring-files)
1. [When to use Jupyter on HPC vs. batch jobs](#jupyter-vs-batch-jobs)

For general information about getting an HPC account, see the [HPC access documentation](../../hpc/access/index.md).

---

## Launching Jupyter Lab on the HPC

1. Go to [https://ondemand-p01.pax.tufts.edu/](https://ondemand-p01.pax.tufts.edu/) and log in with your Tufts SSO credentials
1. Click **Interactive Apps** in the top navigation bar
1. Scroll down to the **Servers** section and select **Jupyter**
1. Fill in the resource request form:
   - **Partition:** Choose `batch` for most work; `gpu` if you need a GPU
   - **Number of hours:** Set an appropriate time limit for your session
   - **Number of CPUs:** one to four is sufficient for most interactive work
   - **Memory (GB):** 8–16 GB is typical; increase for large datasets
   - **Number of GPUs:** Leave at 0 unless you need GPU acceleration
1. Click **Launch**
1. Your request enters the SLURM queue. Once resources are allocated, a **Connect to Jupyter** button appears; click it to open Jupyter Lab in your browser

> **When you are done:** Close your notebooks and return to the **My Interactive Sessions** page on OnDemand. Click **Delete** to end the session and release the resources.

---

## Selecting a Kernel

When you open a new notebook or start a new console, Jupyter asks you to select a kernel. The kernel is the Python environment that runs your code.

By default, several kernels are available corresponding to different Python versions and pre-installed module environments on the HPC. Choose the one that matches your needs.

> **Tip:** If you need specific packages that are not in a default kernel, create your own conda environment and register it as a kernel (see below).

---

## Using a Conda Environment as a Kernel

If you have set up a conda environment on the HPC (see [package management](60-package-management.md)), you can make it available as a Jupyter kernel.

1. **Activate your environment.** Open a terminal on the HPC (via OnDemand: **Clusters > Tufts HPC Shell Access**, or from a VS Code Server session). If conda is not already available, load the module first, then activate your environment:

   ```bash
   module load miniforge/25.3.0
   conda activate myenv
   ```

1. **Install `ipykernel`** in the environment:

   ```bash
   conda install ipykernel
   ```

1. **Register the environment as a kernel:**

   ```bash
   python -m ipykernel install --user --name myenv --display-name "Python (myenv)"
   ```

1. **Restart your Jupyter session** (if it is already running) so the new kernel appears in the kernel selector.

The next time you launch Jupyter on the HPC, your conda environment will be listed as a kernel option under the display name you chose.

### Removing a Kernel

To remove a kernel you no longer need:

```bash
jupyter kernelspec remove myenv
```

---

## Transferring Files

Your HPC home directory is at `/cluster/home/your_username`. Files you create in your Jupyter session are saved there.

### Option 1: Globus (Recommended for Large Files)

[Globus](../../hpc/access/globus/index.md) is the recommended tool for transferring large files to and from the HPC. It is especially useful for datasets that are too large to upload through a browser.

### Option 2: OnDemand File Manager

For small files, use the built-in file manager in OnDemand:

1. Go to **Files > Home Directory** in the OnDemand navigation bar
1. Use the **Upload** button to upload files from your computer
1. Use **Download** to retrieve files

### Option 3: `scp` from the Command Line

For users comfortable with the command line, `scp` (secure copy) is a simple option:

```bash
# Upload a file to the HPC
scp my_data.csv your_username@xfer.pax.tufts.edu:/cluster/home/your_username/

# Download a file from the HPC
scp your_username@xfer.pax.tufts.edu:/cluster/home/your_username/results.csv ./
```

---

## Jupyter vs. Batch Jobs

Jupyter on the HPC is well suited for interactive, exploratory work, but it has limitations:

| Use case                                 | Recommendation  |
| ---------------------------------------- | --------------- |
| Exploratory data analysis                | Jupyter on HPC  |
| Developing and debugging code            | Jupyter on HPC  |
| Short-to-medium analyses (< a few hours) | Jupyter on HPC  |
| Long-running analyses (hours or days)    | SLURM batch job |
| Large-scale runs that can run overnight  | SLURM batch job |
| Automated or scheduled jobs              | SLURM batch job |

Jupyter sessions have time limits and will be terminated if your session time expires. For analyses that might run longer than your session allows, convert your notebook to a `.py` script and submit it as a [SLURM batch job](../../hpc/slurm/index.md) instead.

To convert a notebook to a script:

```bash
jupyter nbconvert --to script my_analysis.ipynb
```

---

## Getting Help

For HPC-related questions, contact Research Technology at **datalab-support@elist.tufts.edu** or submit a request through the [Tufts IT help portal](https://it.tufts.edu/).
