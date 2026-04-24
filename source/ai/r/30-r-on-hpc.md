---
tags: r hpc rstudio
---

# R on the Tufts HPC Cluster

This page provides an overview of how to use R on the Tufts High-Performance Computing (HPC) cluster. The HPC gives you access to much more memory and computing power than a typical laptop or desktop, making it suitable for large datasets, long-running analyses, and parallel computing tasks.

## Getting Access

To use the HPC, you first need an account. Apply at the [Tufts HPC page](https://it.tufts.edu/high-performance-computing). Faculty, staff, and students affiliated with a Tufts research group are eligible.

For general information about using the HPC (including how to log in, transfer files, and submit jobs), see the [HPC documentation](../../hpc/index.md).

## Two Ways to Use R on the HPC

### Option 1: RStudio Server via OnDemand (Recommended for Interactive Analysis)

OnDemand is a web-based interface to the HPC cluster. It includes **RStudio Server**, which gives you a familiar RStudio experience running directly on the cluster, with no command-line setup required.

**To launch RStudio Server:**

1. Go to [https://ondemand-p01.pax.tufts.edu/](https://ondemand-p01.pax.tufts.edu/) and log in with your Tufts SSO credentials
1. Click **Interactive Apps** in the top menu, then select **RStudio Server**
1. Fill in the resource request form:
   - Choose the number of CPU cores and memory appropriate for your analysis
   - Select the R version you need
   - Set a time limit for your session
1. Click **Launch**; your session will start once resources are allocated
1. Click **Connect to RStudio Server** to open the RStudio interface in your browser

> **When you are done:** Exit R properly with `q()`, close the RStudio tab, then return to the **My Interactive Sessions** page on OnDemand and click **Delete** to release the resources. Leaving idle sessions running wastes shared resources.

For more details, including troubleshooting RStudio Server, see the [HPC RStudio guide](../../hpc/application/25-rstudio.md).

### Option 2: R via the Command Line or Batch Jobs

For long-running or automated analyses, you can run R non-interactively by submitting scripts to the SLURM job scheduler. This is more efficient than interactive sessions for jobs that do not require manual input.

#### Loading R

First, connect to the HPC login node:

```bash
ssh your_username@login-prod.pax.tufts.edu
```

Then load the R module:

```bash
module load r/4.4.3
```

Check available versions with:

```bash
module av r
```

#### Running an Interactive R Session

For quick interactive work from the command line, start an interactive compute session:

```bash
srun -p batch -n 2 --mem=4g -t 4:00:00 --pty bash
```

Then start R:

```bash
R
```

For more information on SLURM interactive sessions, see the [SLURM interactive jobs guide](../../hpc/slurm/interactive.md).

#### Submitting an R Batch Job

For long-running analyses, submit your R script as a SLURM batch job. Create a submission script (e.g., `run_analysis.sh`):

```bash
#!/bin/bash
#SBATCH -p batch
#SBATCH -n 4
#SBATCH --mem=16g
#SBATCH -t 12:00:00
#SBATCH -o output_%j.log

module load r/4.4.3
Rscript my_analysis.R
```

Submit the job:

```bash
sbatch run_analysis.sh
```

For more information on SLURM batch jobs, see the [SLURM batch jobs guide](../../hpc/slurm/index.md).

## Managing R Packages on the HPC

Packages you install in your HPC R sessions are stored in your home directory and are separate from your local R installation. Install packages the same way you would locally, from an R or RStudio session on the HPC:

```r
install.packages("tidyverse")
```

For reproducible package environments, you can use **renv** on the HPC just as you would locally. See the [RStudio Desktop setup guide](20-rstudio-desktop.md) for an introduction to renv.

> **Note:** Some packages require compiled dependencies that are available on the HPC as modules. If you encounter installation errors, contact Research Technology support at datalab-support@elist.tufts.edu.

## Data Storage on the HPC

Files on the HPC are stored in your home directory (`/cluster/home/your_username`) and in research storage allocations. For information about data storage options and policies, see the [HPC storage documentation](../../hpc/access/index.md).

> **Data security:** Check with your IRB or Data Use Agreement before transferring sensitive data to the HPC. For guidance, consult the [Tufts data storage finder](https://access.tufts.edu/data-finder).

## Getting Help

For HPC-related questions, contact Research Technology at **datalab-support@elist.tufts.edu** or submit a request through the [Tufts IT help portal](https://it.tufts.edu/).
