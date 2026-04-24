---
tags: matlab hpc slurm batch ondemand
---

# MATLAB on the Tufts HPC Cluster

The Tufts HPC cluster provides two ways to run MATLAB:

- **MATLAB Server via OnDemand** (browser-based): the recommended interactive approach — no X11 or VDI required
- **Batch mode via Slurm**: the recommended approach for long or computationally intensive analyses

```{attention}

**Do not run MATLAB on the login nodes.** Login nodes are shared by all users for file management
and job submission only. All MATLAB work — batch or interactive — must use compute nodes
via Slurm or the OnDemand MATLAB Server.
```

---

## Prerequisites

### 1. Cluster access

You must have an active Tufts HPC account. If you do not have one, request access at the HPC Account Request page.

### 2. Tufts MATLAB license

All Tufts faculty, staff, and students are covered by the campus license. See [access.tufts.edu/matlab](https://access.tufts.edu/matlab) for details.

### 3. VPN (off-campus users)

If you are connecting from off campus, connect to the [Tufts VPN](https://access.tufts.edu/vpn) before proceeding.

### 4. Available MATLAB versions

The following versions are available on the cluster:

```
matlab/R2025a
matlab/2023b
matlab/2023b-beta
matlab/2023a
matlab/2022a
matlab/2021a
matlab/2020b
matlab/2020a
```

Load a specific version:

```bash
module load matlab/R2025a
```

To check what is available or get details on a version:

```bash
module spider matlab
module spider matlab/R2025a
```

---

## Method 1: MATLAB Server via OnDemand (Recommended for Interactive Work)

The OnDemand MATLAB Server runs a full MATLAB session on a compute node and streams the interface to your browser — no X11 forwarding or local software required.

**Step 1 — Log in to OnDemand:**

Navigate to [https://ondemand-prod.pax.tufts.edu](https://ondemand-prod.pax.tufts.edu) and log in with your Tufts credentials.

**Step 2 — Launch MATLAB Server:**

Select **MATLAB Server** from the **Interactive Apps** menu.

**Step 3 — Configure your session:**

Fill out the launch form with your resource requirements:

| Field           | Guidance                                                      |
| --------------- | ------------------------------------------------------------- |
| MATLAB version  | Select from available versions (e.g. `matlab/R2025a`)         |
| Number of hours | Estimate your session length; session ends when time expires  |
| Number of cores | 1 for most work; increase for Parallel Computing Toolbox jobs |
| Memory (GB)     | 8–16 GB for most analyses; increase for large datasets        |
| Partition       | `interactive` for short sessions; `batch` for longer work     |

Click **Launch**. Your session will be queued and start when resources are available.

**Step 4 — Connect:**

Once the session is ready, click **Connect to MATLAB Server**. The full MATLAB interface opens in your browser.

**Step 5 — End your session:**

When finished, exit MATLAB and click **Delete** on the session card in OnDemand to release the compute node for other users.

---

## Method 2: Batch Mode via Slurm (Recommended for Long Jobs)

Batch mode runs MATLAB non-interactively using a `.m` script file. This is the preferred approach for analyses that take more than a few minutes, require large memory, or need to run unattended.

### Step 1 — Write your MATLAB script

Create your MATLAB script on the cluster. For example, `analysis.m`:

```matlab
% analysis.m — example MATLAB script
data = load('/cluster/tufts/mylab/myutln/data/mydata.mat');

results = mean(data.values);
disp(results);

save('/cluster/tufts/mylab/myutln/output/results.mat', 'results');
```

### Step 2 — Write a Slurm batch script

Create a Slurm submission script, e.g. `run_matlab.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=matlab_analysis
#SBATCH --output=matlab_%j.out
#SBATCH --error=matlab_%j.err
#SBATCH --time=04:00:00          # adjust to your expected run time
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1        # increase for parfor or Parallel Computing Toolbox
#SBATCH --mem=16G                # adjust to your dataset size
#SBATCH --partition=batch

module purge
module load matlab/R2025a

# -nodisplay -nodesktop -nosplash: suppress GUI
# -nojvm: skip Java VM (faster startup; omit if you need figure/plot output)
matlab -nodisplay -nodesktop -nosplash -nojvm -r "run('analysis.m'); exit"
```

```{note}
Omit `-nojvm` if your script generates figures or uses toolboxes that require the Java VM.
Use `-batch` instead of `-r` in newer MATLAB versions (R2019b+) — it handles errors and exit codes more reliably:

`matlab -batch "run('analysis.m')"`
```

### Step 3 — Submit the job

```bash
sbatch run_matlab.sh
```

### Step 4 — Check output

Monitor job status:

```bash
squeue -u $USER
```

Review output once the job completes:

```bash
cat matlab_<jobid>.out
cat matlab_<jobid>.err
```

---

### Running Multiple Scripts as a Job Array

If you have several independent MATLAB scripts to run (e.g. one per dataset), use a Slurm job array:

```bash
#!/bin/bash
#SBATCH --job-name=matlab_array
#SBATCH --output=matlab_%A_%a.out
#SBATCH --error=matlab_%A_%a.err
#SBATCH --array=1-10              # one task per script
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=batch

module purge
module load matlab/R2025a

matlab -batch "run('job_${SLURM_ARRAY_TASK_ID}.m')"
```

---

### Using the Parallel Computing Toolbox

MATLAB's Parallel Computing Toolbox (`parfor`, `parfeval`) can use multiple cores on a single node. Request the cores you need in your Slurm script and configure a local parallel pool in your MATLAB code:

```bash
#SBATCH --cpus-per-task=8
```

```matlab
% In your .m script — create a pool matching the allocated cores
parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));

parfor i = 1:100
    % your parallel work here
end

delete(gcp('nocreate'));
```

---

## Quick Reference

| Task                       | Command                                                             |
| -------------------------- | ------------------------------------------------------------------- |
| Check available versions   | `module spider matlab`                                              |
| Load MATLAB                | `module load matlab/R2025a`                                         |
| Run script in batch mode   | `matlab -batch "run('script.m')"`                                   |
| Run script (older syntax)  | `matlab -nodisplay -nodesktop -nosplash -r "run('script.m'); exit"` |
| Submit batch job           | `sbatch run_matlab.sh`                                              |
| Monitor job                | `squeue -u $USER`                                                   |
| Launch interactive session | OnDemand → Interactive Apps → MATLAB Server                         |

---

## Getting Help

- **Tufts MATLAB license:** [access.tufts.edu/matlab](https://access.tufts.edu/matlab)
- **Email:** [tts-research@tufts.edu](mailto:tts-research@tufts.edu)
- **MATLAB documentation:** [mathworks.com/help/matlab](https://www.mathworks.com/help/matlab/)
- **MATLAB Parallel Computing Toolbox:** [mathworks.com/help/parallel-computing](https://www.mathworks.com/help/parallel-computing/)

> Linked external resources are not affiliated with or endorsed by Tufts University.
