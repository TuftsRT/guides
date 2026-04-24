---
tags: sas data-science statistics hpc
---

# SAS on the Tufts HPC Cluster

SAS (Statistical Analysis System) is a comprehensive software suite for advanced analytics, data management, and statistical analysis. This guide covers how to run SAS on the Tufts HPC cluster in two modes:

- **Batch mode** (command line): the recommended approach for long or computationally intensive analyses
- **Interactive GUI mode** (SAS Display Manager via X11 forwarding): useful for exploration, visualization, and debugging

```{attention}

**Do not run SAS on the login nodes.** Login nodes are shared by all users for file management
and job submission only. All SAS work — batch or interactive — must be submitted to compute nodes
via Slurm.
```

---

## Prerequisites

### 1. Cluster access

You must have an active Tufts HPC account. If you do not have one, request access at the HPC Account Request page.

### 2. VPN (off-campus users)

If you are connecting to the Tufts HPC from off campus, connect to the [Tufts VPN](https://access.tufts.edu/vpn) before proceeding.

### 3. Load SAS

SAS is available via the deprecated module tree. Login via the [Open OnDemand portal](https://ondemand-p01.pax.tufts.edu), navigate to Clusters > Tufts HPC Shell Access, and load SAS with:

```bash
module load modtree/deprecated SAS/9.4
```

---

## Method 1: Batch Mode (Recommended)

Batch mode runs SAS non-interactively using a **SAS program file** (`.sas`). This is the preferred method for any analysis that takes more than a few minutes, and is the only reliable approach for long or memory-intensive jobs.

### Step 1 — Write your SAS program

Create your SAS script locally or on the cluster. For example, `analysis.sas`:

```sas
/* analysis.sas — example SAS program */
libname mydata "/cluster/tufts/mylab/myutln/data";

proc contents data=mydata.myfile;
run;

proc means data=mydata.myfile;
run;

proc reg data=mydata.myfile;
  model y = x1 x2 x3;
run;
```

### Step 2 — Write a Slurm batch script

Create a Slurm submission script, e.g. `run_sas.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=sas_analysis
#SBATCH --output=sas_%j.out
#SBATCH --error=sas_%j.err
#SBATCH --time=04:00:00          # adjust to your expected run time
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1        # increase for SAS/STAT parallel procedures
#SBATCH --mem=16G                # adjust to your dataset size
#SBATCH --partition=batch

# Load SAS
module load modtree/deprecated SAS/9.4

# Run SAS in batch mode
# -noterminal suppresses the Display Manager; output goes to analysis.log and analysis.lst
sas analysis.sas -noterminal
```

### Step 3 — Submit the job

```bash
sbatch run_sas.sh
```

### Step 4 — Check output

SAS batch mode produces two output files:

- `analysis.log` — the SAS log (warnings, errors, timing)
- `analysis.lst` — the procedure output (tables, statistics)

Review them once the job completes:

```bash
cat analysis.log
cat analysis.lst
```

You can also monitor job status with:

```bash
squeue -u $USER
```

---

### Running Multiple SAS Programs as a Job Array

If you have several independent SAS programs to run (e.g. processing one dataset per file), use a Slurm job array rather than submitting individual jobs:

```bash
#!/bin/bash
#SBATCH --job-name=sas_array
#SBATCH --output=sas_%A_%a.out
#SBATCH --error=sas_%A_%a.err
#SBATCH --array=1-10              # one task per SAS program
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=batch

module load modtree/deprecated SAS/9.4

# Each task runs the SAS program matching its array index
sas "job_${SLURM_ARRAY_TASK_ID}.sas" -noterminal
```

---

## Method 2: Interactive GUI Mode (SAS Display Manager via X11)

The SAS Display Manager launches the full SAS graphical interface on a compute node and forwards the display back to your local screen via X11. This is useful for:

- Exploratory data analysis
- Building and testing SAS programs interactively
- Viewing output and graphs in real time

```{attention}
Interactive SAS GUI sessions consume cluster resources for their entire duration.
Please exit SAS and end your session when you are done. Do not leave idle SAS sessions
running overnight.
```

### X11 Setup by Operating System

Before launching the SAS GUI, your local machine must be able to receive X11 graphics. See the X11 forwarding setup guide here: https://rtguides.it.tufts.edu/hpc/access/20-cli.html

---

### Launching the SAS Display Manager

Once connected with X11 forwarding, **do not run SAS on the login node.** Request an interactive compute session first, then load SAS and launch the GUI.

**Step 1 — Request an interactive compute node:**

```bash
srun --pty --partition=interactive --time=02:00:00 --mem=8G --cpus-per-task=1 bash
```

Adjust `--time` and `--mem` to match your needs. When your session is granted, your prompt will change to reflect the compute node (e.g. `[your_utln@c1234 ~]$`).

**Step 2 — Load SAS and launch the GUI:**

```bash
module load modtree/deprecated SAS/9.4
sas &
```

The `&` runs SAS in the background so you retain use of the terminal. The SAS Display Manager will appear on your local screen via X11.

**Step 3 — Work interactively, then exit cleanly:**

When you are finished, exit SAS via the GUI (`File > Exit`) or submit `endsas;` in the Program Editor window. Then type `exit` in the terminal to release the compute node.

---

### Troubleshooting X11 / SAS GUI

**"cannot open display" or blank window**

- Confirm you used `-X` or `-Y` when SSH-ing to the cluster.
- On macOS, verify XQuartz is running (`Applications > Utilities > XQuartz`).
- Check that your `DISPLAY` variable is set: `echo $DISPLAY` — it should return something like `localhost:10.0`. If it is empty, your SSH connection does not have X11 forwarding active.

**GUI is very slow or laggy**

- X11 over SSH can be slow on high-latency connections. Use `-Y` (trusted forwarding) instead of `-X`, or consider using the OnDemand HPC Desktop as an alternative — it runs a full Linux desktop in your browser without requiring X11 setup.

**SAS crashes or fails to start**

- Try running `sas -noterminal` first to confirm the module loaded correctly.
- Ensure you are on a compute node, not a login node.

---

## Alternative: OnDemand HPC Desktop

If X11 forwarding is proving difficult to configure — particularly on Windows or over VPN — the OnDemand HPC Desktop provides a full Linux desktop environment in your browser. From the desktop, you can open a terminal, load SAS, and run the GUI without any local X11 software.

Access OnDemand at: https://ondemand-prod.pax.tufts.edu

Navigate to **Interactive Apps → HPC Desktop**, request a session, and once it launches, open a terminal within the desktop and follow the same steps above (load module, run `sas &`).

---

## Installing SAS Macros and Custom Formats

User-written SAS macros and custom formats can be stored in your home or project directory and referenced in your SAS programs using `%include` or `libname`/`options sasautos=`.

To make macros available to all your SAS jobs, add the following to your SAS program (or a site-level `autoexec.sas`):

```sas
options sasautos=("/cluster/tufts/mylab/myutln/sas_macros" sasautos);
```

Custom format catalogs should be stored in a dedicated library and referenced with `libname` before any procedure that uses them.

---

## Quick Reference

| Task                             | Command                                  |
| -------------------------------- | ---------------------------------------- |
| Check available versions         | `module spider sas`                      |
| Load SAS                         | `module load modtree/deprecated SAS/9.4` |
| Run SAS program in batch mode    | `sas analysis.sas -noterminal`           |
| Launch SAS GUI (after X11 login) | `sas &`                                  |
| View SAS log                     | `cat analysis.log`                       |
| View SAS output listing          | `cat analysis.lst`                       |
| Monitor job status               | `squeue -u $USER`                        |

---

## Getting Help

For questions about running SAS on the Tufts HPC cluster:

- **Email:** [tts-research@tufts.edu](mailto:tts-research@tufts.edu)
- **SAS documentation:** [documentation.sas.com](https://documentation.sas.com/)
- **SAS communities:** [communities.sas.com](https://communities.sas.com/)

> Linked external resources are not affiliated with or endorsed by Tufts University.
