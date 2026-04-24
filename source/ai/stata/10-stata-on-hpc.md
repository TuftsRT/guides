# Stata on the Tufts HPC Cluster

Stata is a complete, integrated statistical software package for data analysis, data management, and graphics. This guide covers how to run Stata on the Tufts HPC cluster in two modes:

- **Batch mode** (command line): the recommended approach for long or computationally intensive analyses
- **Interactive GUI mode** (xstata via X11 forwarding): useful for exploration, visualization, and debugging

We will also discuss topics such as:

- **Parallel processing**: using Stata/MP or the `parallel` package for multi-core jobs
- **Large dataset strategies**: compression, subsetting, chunking, and memory-efficient packages

```{attention}

**Do not run Stata on the login nodes.** Login nodes are shared by all users for file management
and job submission only. All Stata work — batch or interactive — must be submitted to compute nodes
via Slurm.
```

---

## Prerequisites

### 1. Cluster access

You must have an active Tufts HPC account. If you do not have one, request access at the HPC Account Request page.

### 2. VPN (off-campus users)

If you are connecting to the Tufts HPC from off campus, connect to the [Tufts VPN](https://access.tufts.edu/vpn) before proceeding.

### 3. Check available Stata versions

Login via the [Open On Demand portal](https://ondemand-p01.pax.tufts.edu), and navigate to the Clusters > Tufts HPC Shell Access.

Once logged in, check which versions are installed:

```bash
module spider stata
```

Load the default (or a specific) version:

```bash
module load stata
# or, for a specific version:
module load stata/19
```

---

## Method 1: Batch Mode (Recommended)

Batch mode runs Stata non-interactively using a **do-file** (a plain text file of Stata commands). This is the preferred method for any analysis that takes more than a few minutes, and is the only reliable approach for long or memory-intensive jobs.

### Step 1 — Write your do-file

Create your Stata script locally or on the cluster. For example, `analysis.do`:

```stata
* analysis.do — example Stata do-file
use "/cluster/tufts/mylab/myutln/data/mydata.dta", clear

summarize
regress y x1 x2 x3

log using results.log, replace
estimates table, b se
log close
```

### Step 2 — Write a Slurm batch script

Create a Slurm submission script, e.g. `run_stata.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=stata_analysis
#SBATCH --output=stata_%j.out
#SBATCH --error=stata_%j.err
#SBATCH --time=04:00:00          # adjust to your expected run time
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1        # use 1 CPU for Stata/SE; see Stata/MP note below
#SBATCH --mem=16G                # adjust to your dataset size
#SBATCH --partition=batch

# Load Stata
module load stata

# Run Stata in batch (non-interactive) mode
# -b flag suppresses the GUI and writes output to analysis.log
stata -b do analysis.do
```

### Step 3 — Submit the job

```bash
sbatch run_stata.sh
```

### Step 4 — Check output

Stata batch mode automatically writes output to a log file named after your do-file (`analysis.log`). Review it once the job completes:

```bash
cat analysis.log
```

You can also monitor job status with:

```bash
squeue -u $USER
```

---

### Running Multiple Do-Files as a Job Array

If you have several independent do-files to run (e.g. processing one dataset per file), use a Slurm job array rather than submitting individual jobs:

```bash
#!/bin/bash
#SBATCH --job-name=stata_array
#SBATCH --output=stata_%A_%a.out
#SBATCH --error=stata_%A_%a.err
#SBATCH --array=1-10              # one task per do-file
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=batch

module load stata

# Each task runs the do-file matching its array index
stata -b do "job_${SLURM_ARRAY_TASK_ID}.do"
```

---

## Parallel Processing in Stata

### Stata/MP

Stata/MP is the multi-processor edition of Stata. Many core estimation and data management commands are parallelized automatically — you do not need to modify your do-file. To use Stata/MP, request multiple CPUs in your Slurm script and invoke `stata-mp` instead of `stata`:

```bash
#!/bin/bash
#SBATCH --job-name=stata_mp
#SBATCH --output=stata_%j.out
#SBATCH --error=stata_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4        # request cores for Stata/MP
#SBATCH --mem=32G
#SBATCH --partition=batch

module load stata

# Tell Stata/MP how many cores to use
stata-mp -b -n 4 do analysis.do
```

The `-n` flag sets the number of processors. As a rule of thumb, start with 4 cores — Stata/MP scales well up to 4 cores for most procedures, but returns diminish beyond that for many analyses. Check [Stata's MP benchmark data](https://www.stata.com/statamp/) to see whether your specific commands benefit from more cores before requesting them.

```{note}
Not all Stata commands are parallelized in Stata/MP. Commands that are not MP-aware will run on a single core regardless of how many you request. Requesting more cores than your analysis can use wastes your fairshare allocation.
```

### The `parallel` Package

The community-contributed [`parallel`](https://github.com/gvegayon/parallel) package (by George Vega Yon) enables embarrassingly parallel workflows by spawning multiple Stata instances and splitting work across them. This is useful for tasks like bootstrap resampling, simulation, or applying the same analysis to many independent subgroups.

**Install once** (from an interactive or login-node session):

```stata
. ssc install parallel, replace
. parallel initialize 4   // set default number of cores
```

**Example — parallel bootstrap:**

```stata
parallel initialize 4

* Run 1000 bootstrap replications across 4 Stata instances
parallel bs, reps(1000): regress y x1 x2 x3
```

**Example — apply analysis to independent groups:**

```stata
parallel initialize 4

* Runs the block in parallel for each value of groupvar
parallel, by(groupvar): {
    regress y x1 x2 x3
    estimates save "results_`groupvar'", replace
}
```

In a Slurm batch script, match `--cpus-per-task` to the number of cores you pass to `parallel initialize`:

```bash
#SBATCH --cpus-per-task=4
...
stata -b do parallel_analysis.do
```

---

## Working with Large Datasets

Stata loads data entirely into RAM by default. When datasets approach or exceed available memory, use the strategies below to keep jobs running efficiently.

### Check and reduce memory usage

Always run `compress` before a memory-intensive analysis. It automatically downcasts variables to the smallest storage type that preserves their values — often cutting dataset size by 30–60%:

```stata
use "/cluster/tufts/mylab/myutln/data/mydata.dta", clear
compress
save "/cluster/tufts/mylab/myutln/data/mydata_compressed.dta", replace
```

Check memory consumption with:

```stata
describe, short    // shows dataset size in memory
memory             // shows Stata's current memory allocation
```

### Load only what you need

Avoid loading variables or observations you do not use:

```stata
* Load a subset of variables
use id year income wage using "mydata.dta", clear

* Load a subset of observations matching a condition
use "mydata.dta" if year >= 2010, clear

* Load a range of observations (useful for testing on a slice)
use in 1/100000 using "mydata.dta", clear
```

Combining these:

```stata
use id year income if year >= 2010 using "mydata.dta", clear
```

### Process data in chunks

For datasets too large to fit in memory even after compression, split the work into chunks, process each chunk separately, and combine results:

```stata
* chunk_analysis.do
* Called with: stata -b do chunk_analysis.do <start> <end>

local start = `1'
local end   = `2'

use in `start'/`end' using "/cluster/tufts/mylab/myutln/data/bigdata.dta", clear
compress

* ... your analysis ...

save "/cluster/tufts/mylab/myutln/output/results_`start'_`end'.dta", replace
```

Submit as a job array, one task per chunk:

```bash
#!/bin/bash
#SBATCH --array=0-9
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=batch

module load stata

# Each task processes 1,000,000 rows; adjust to your dataset size
START=$(( SLURM_ARRAY_TASK_ID * 1000000 + 1 ))
END=$(( (SLURM_ARRAY_TASK_ID + 1) * 1000000 ))

stata -b do chunk_analysis.do $START $END
```

Then combine the chunk results in a final do-file:

```stata
* combine_results.do
clear
forvalues i = 0/9 {
    local start = `i' * 1000000 + 1
    local end   = (`i' + 1) * 1000000
    append using "/cluster/tufts/mylab/myutln/output/results_`start'_`end'.dta"
}
save "/cluster/tufts/mylab/myutln/output/combined_results.dta", replace
```

### Use memory-efficient community packages

The [`gtools`](https://gtools.readthedocs.io) package provides drop-in replacements for many slow, memory-heavy Stata commands (`collapse`, `egen`, `xtile`, `levelsof`, etc.) using optimized C plugins. It is significantly faster and uses less peak memory than the built-in equivalents:

```stata
ssc install gtools, replace

* Instead of: collapse (mean) income wage, by(state year)
gcollapse (mean) income wage, by(state year)

* Instead of: egen group_mean = mean(income), by(state)
gegen group_mean = mean(income), by(state)
```

For panel and time-series data with many groups, [`ftools`](https://github.com/sergiocorreia/ftools) offers similar speedups for `fcollapse`, `fsort`, and `fmerge`.

---

## Method 2: Interactive GUI Mode (xstata via X11)

`xstata` launches the full Stata graphical interface on a compute node and forwards the display back to your local screen via X11. This is useful for:

- Exploratory data analysis
- Building and testing do-files interactively
- Viewing graphs and results in real time

```{attention}
xstata sessions are **interactive jobs** and consume cluster resources for their entire duration.
Please exit Stata and end your session when you are done. Do not leave idle xstata sessions
running overnight.
```

### X11 Setup by Operating System

Before launching xstata, your local machine must be able to receive X11 graphics. Check out the details on x11 forwarding here: https://rtguides.it.tufts.edu/hpc/access/20-cli.html

---

### Launching xstata

Once connected with X11 forwarding, **do not run xstata on the login node.** Request an interactive compute session first, then load Stata and launch the GUI.

**Step 1 — Request an interactive compute node:**

```bash
srun --pty --partition=interactive --time=02:00:00 --mem=8G --cpus-per-task=1 bash
```

Adjust `--time` and `--mem` to match your needs. When your session is granted, your prompt will change to reflect the compute node (e.g. `[your_utln@c1234 ~]$`).

**Step 2 — Load Stata and launch the GUI:**

```bash
module load stata
xstata &
```

The `&` runs xstata in the background so you retain use of the terminal. The Stata graphical interface will appear on your local screen via X11.

**Step 3 — Work interactively, then exit cleanly:**

When you are finished, exit Stata via the GUI (`File > Exit`) or type `exit, clear` in the Stata command window. Then type `exit` in the terminal to release the compute node.

---

### Troubleshooting X11 / xstata

**"cannot open display" or blank window**

- Confirm you used `-X` or `-Y` when SSH-ing to the cluster.
- On macOS, verify XQuartz is running (`Applications > Utilities > XQuartz`).
- Check that your `DISPLAY` variable is set: `echo $DISPLAY` — it should return something like `localhost:10.0`. If it is empty, your SSH connection does not have X11 forwarding active.

**GUI is very slow or laggy**

- X11 over SSH can be slow on high-latency connections. Use `-Y` (trusted forwarding) instead of `-X`, or consider using the OnDemand HPC Desktop as an alternative — it runs a full Linux desktop in your browser without requiring X11 setup.

**xstata crashes immediately**

- Try running `stata` (the command-line version) first to confirm the module loaded correctly.
- Ensure you are on a compute node, not a login node.

---

## Alternative: OnDemand HPC Desktop

If X11 forwarding is proving difficult to configure — particularly on Windows or over VPN — the OnDemand HPC Desktop provides a full Linux desktop environment in your browser. From the desktop, you can open a terminal, load Stata, and run `xstata` without any local X11 software.

Access OnDemand at: https://ondemand-prod.pax.tufts.edu

Navigate to **Interactive Apps → HPC Desktop**, request a session, and once it launches, open a terminal within the desktop and follow the same steps as above (load module, run `xstata`).

---

## Installing Stata Packages (ado files)

Stata community-contributed packages (installed via `ssc install` or `net install`) are written to your personal ado directory, typically `~/ado/plus/`. These persist across sessions and do not need to be reinstalled each time.

Install packages from an **interactive session on a login node** or a compute node with internet access:

```stata
. ssc install estout, replace
. ssc install outreg2, replace
```

If compute nodes do not have internet access, install packages from the login node (without Slurm) during a quick interactive terminal session:

```bash
# On the login node only — lightweight package installation, not analysis
module load stata
stata
```

Then run `ssc install` within the Stata prompt. Installed packages will be available to all subsequent batch and interactive jobs.

---

## Quick Reference

| Task                                    | Command                             |
| --------------------------------------- | ----------------------------------- |
| Check available versions                | `module spider stata`               |
| Load Stata                              | `module load stata`                 |
| Run do-file in batch mode               | `stata -b do yourscript.do`         |
| Run Stata/MP with N cores               | `stata-mp -b -n 4 do yourscript.do` |
| Launch Stata GUI (after X11 login)      | `xstata &`                          |
| Launch Stata command line interactively | `stata`                             |
| Compress dataset in memory              | `compress`                          |
| Check dataset memory usage              | `describe, short`                   |
| Load subset of variables                | `use varlist using file.dta, clear` |
| Load subset of observations             | `use file.dta if condition, clear`  |
| View batch log output                   | `cat yourscript.log`                |
| Install an ado package                  | `ssc install <package>, replace`    |

---

## Getting Help

For questions about running Stata on the Tufts HPC cluster:

- **Email:** [tts-research@tufts.edu](mailto:tts-research@tufts.edu)
- **Stata documentation:** [stata.com/features/](https://www.stata.com/features/)
- **Stata user community:** [statalist.org](https://www.statalist.org/)

> Linked external resources are not affiliated with or endorsed by Tufts University.
