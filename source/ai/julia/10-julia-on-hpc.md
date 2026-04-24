---
tags: julia hpc slurm batch ondemand jupyter
---

# Julia on the Tufts HPC Cluster

Julia is a high-level, high-performance programming language designed for numerical and scientific computing. It combines the ease of use of Python with speeds approaching C, making it well suited for simulations, optimization, differential equations, machine learning, and large-scale data analysis.

This guide covers:

- **Batch mode via Slurm**: the recommended approach for long or computationally intensive analyses
- **Interactive sessions via OnDemand**: for exploration, prototyping, and Jupyter notebooks
- **Package management**: Julia's built-in package manager for reproducible environments

```{attention}

**Do not run Julia on the login nodes.** Login nodes are shared by all users for file management
and job submission only. All Julia work — batch or interactive — must use compute nodes
via Slurm or OnDemand.
```

---

## Prerequisites

### 1. Cluster access

You must have an active Tufts HPC account. If you do not have one, request access at the HPC Account Request page.

### 2. VPN (off-campus users)

If you are connecting from off campus, connect to the [Tufts VPN](https://access.tufts.edu/vpn) before proceeding.

### 3. Load Julia

The following Julia versions are available on the cluster:

```
julia/1.11.0   ← recommended
julia/1.7.1
julia/1.6.0
julia/1.5.3
julia/0.3
```

Load the recommended version:

```bash
module load julia/1.11.0
```

To get details on a specific version or check for new additions:

```bash
module spider julia
module spider julia/1.11.0
```

---

## Method 1: Batch Mode via Slurm (Recommended for Long Jobs)

Batch mode runs Julia non-interactively using a `.jl` script. This is the preferred approach for analyses that take more than a few minutes or require significant memory.

### Step 1 — Write your Julia script

Create your Julia script on the cluster. For example, `analysis.jl`:

```julia
# analysis.jl — example Julia script
using Statistics, DelimitedFiles

data = readdlm("/cluster/tufts/mylab/myutln/data/mydata.csv", ',', Float64, '\n')

println("Mean: ", mean(data))
println("Std:  ", std(data))

writedlm("/cluster/tufts/mylab/myutln/output/results.csv", data, ',')
```

### Step 2 — Write a Slurm batch script

Create a Slurm submission script, e.g. `run_julia.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=julia_analysis
#SBATCH --output=julia_%j.out
#SBATCH --error=julia_%j.err
#SBATCH --time=04:00:00          # adjust to your expected run time
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1        # increase for multi-threaded Julia
#SBATCH --mem=16G                # adjust to your dataset size
#SBATCH --partition=batch

module load julia/1.11.0

julia analysis.jl
```

### Step 3 — Submit the job

```bash
sbatch run_julia.sh
```

### Step 4 — Check output

```bash
squeue -u $USER          # monitor job status
cat julia_<jobid>.out    # review output
cat julia_<jobid>.err    # review errors
```

---

### Multi-threaded Julia

Julia has built-in multi-threading via `Threads.@threads`. To use it, request multiple CPUs in Slurm and set the `JULIA_NUM_THREADS` environment variable:

```bash
#SBATCH --cpus-per-task=8

module load julia/1.11.0
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
julia analysis.jl
```

In your Julia script:

```julia
println("Running with $(Threads.nthreads()) threads")

Threads.@threads for i in 1:100
    # parallel work here
end
```

---

### Running Multiple Scripts as a Job Array

For independent jobs (e.g. one per dataset), use a Slurm job array:

```bash
#!/bin/bash
#SBATCH --job-name=julia_array
#SBATCH --output=julia_%A_%a.out
#SBATCH --error=julia_%A_%a.err
#SBATCH --array=1-10
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=batch

module load julia/1.11.0

julia analysis.jl $SLURM_ARRAY_TASK_ID
```

Then in your Julia script, read the task index from the command line:

```julia
task_id = parse(Int, ARGS[1])
# use task_id to select your dataset or parameter set
```

---

## Method 2: Interactive Session via OnDemand

For exploration and prototyping, request an interactive compute session through the OnDemand shell or Jupyter interface.

### Julia in a terminal

1. Log in to [OnDemand](https://ondemand-prod.pax.tufts.edu) and open a shell on a compute node via **Clusters > Tufts HPC Shell Access**.
1. Request an interactive session:

```bash
srun --pty --partition=interactive --time=02:00:00 --mem=8G --cpus-per-task=1 bash
```

3. Load Julia and launch the REPL:

```bash
module load julia/1.11.0
julia
```

### Julia in Jupyter

You can run Julia in a Jupyter notebook on the HPC cluster by registering a Julia kernel.

**Step 1 — Install IJulia (one-time setup):**

From an interactive session on a compute node:

```bash
module load julia/1.11.0
julia
```

Then inside the Julia REPL:

```julia
using Pkg
Pkg.add("IJulia")
using IJulia
installkernel("Julia")
```

**Step 2 — Launch Jupyter via OnDemand:**

Navigate to **Interactive Apps → Jupyter** in OnDemand, start a session, and select the created **Julia** kernel when creating a new notebook.

---

## Package Management

Julia has a built-in package manager (`Pkg`) that creates per-project environments, similar to Python's conda environments or R's renv. This makes Julia projects reproducible and self-contained.

### Adding packages

From the Julia REPL, enter package mode with `]`:

```julia
julia> ]
(@v1.10) pkg> add DataFrames CSV Statistics Plots
```

Or from a script:

```julia
using Pkg
Pkg.add(["DataFrames", "CSV", "Statistics", "Plots"])
```

### Project environments

For reproducible research, create a project-specific environment:

```bash
mkdir /cluster/tufts/mylab/myutln/myproject
cd /cluster/tufts/mylab/myutln/myproject
julia --project=.
```

Inside the REPL:

```julia
]
(myproject) pkg> add DataFrames CSV
```

This creates `Project.toml` and `Manifest.toml` files that pin all dependencies. To restore the environment on another machine or job:

```julia
using Pkg
Pkg.instantiate()
```

In a Slurm script, activate the project environment at startup:

```bash
julia --project=/cluster/tufts/mylab/myutln/myproject analysis.jl
```

### Package storage

By default, Julia installs packages to `~/.julia/`. On the cluster, this is your home directory. If you install many large packages, monitor your home quota:

```bash
quota -s
```

## You can also store these files in your lab storage directory, if you have one. We recommend this.

## Quick Reference

| Task                              | Command                               |
| --------------------------------- | ------------------------------------- |
| Check available versions          | `module spider julia`                 |
| Load Julia                        | `module load julia/1.11.0`            |
| Run a script                      | `julia script.jl`                     |
| Run with multiple threads         | `JULIA_NUM_THREADS=8 julia script.jl` |
| Start interactive REPL            | `julia`                               |
| Add a package                     | `] add PackageName`                   |
| Activate a project environment    | `julia --project=/path/to/project`    |
| Instantiate (restore) environment | `Pkg.instantiate()`                   |
| Submit batch job                  | `sbatch run_julia.sh`                 |
| Monitor job                       | `squeue -u $USER`                     |

---

## Getting Help

- **Email:** [tts-research@tufts.edu](mailto:tts-research@tufts.edu)
- **Julia documentation:** [docs.julialang.org](https://docs.julialang.org/)
- **Julia Discourse (community):** [discourse.julialang.org](https://discourse.julialang.org/)
- **Julia package registry:** [juliahub.com](https://juliahub.com/)

> Linked external resources are not affiliated with or endorsed by Tufts University.
