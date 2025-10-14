# How Many Cores Do You Need for R Jobs on HPC?

Efficiently using multiple cores can dramatically reduce runtime for R jobs on the Tufts HPC cluster. This guide explains how to estimate the number of cores you need, detect available cores, and enable parallel computing in R.

---

## 1. How to Know How Many Cores You Need

There’s no single “correct” number of cores — it depends on the type of task:

| Task Type                                                             | Typical Core Need     | Notes                                                |
| --------------------------------------------------------------------- | --------------------- | ---------------------------------------------------- |
| Basic data manipulation, plotting, regression                         | 1 core                | R is single-threaded by default                      |
| Linear algebra (matrix ops, e.g., `lm()`, `prcomp()`, `svd()`)        | 1–8 cores             | Some functions use multi-threaded BLAS automatically |
| Parallel processing (`parallel`, `future`, `foreach`, `BiocParallel`) | N cores (you specify) | You control how many cores via parameters            |
| Big data tools (`data.table`, `Seurat`, `DESeq2`, `edgeR`)            | 1–8 cores             | Some have parallel options you can enable manually   |

**Rule of thumb for HPC jobs**

- Start with **4–8 cores** unless you know it scales well.
- Monitor runtime and CPU usage (`sacct`, `squeue`, or `top`), then adjust.

---

## 2. Check How Many Cores Are Available

```r
parallel::detectCores()
```

This tells you the total number of logical cores visible to R.

On a cluster, R will *see* all cores on the node — but you should only use the ones allocated to your job.

---

## 3. Built-in Multi-Core Functionality in R

R can take advantage of multiple cores in two main ways:

### A. Multi-Threaded Math Libraries (BLAS)

If your R installation is linked against **OpenBLAS**, **MKL**, or **BLIS**, many matrix operations automatically use multiple cores.

Check your current BLAS library:

```
sessionInfo()
```

You’ll see something like:

```
linked to BLAS: /usr/lib/x86_64-linux-gnu/openblas/libblas.so
```

Control the number of threads:

```
library(RhpcBLASctl)
blas_get_num_procs()
blas_set_num_threads(4)
```

---

### B. R’s Parallel Frameworks

| Package                    | Typical Function            | Example                                                    |
| -------------------------- | --------------------------- | ---------------------------------------------------------- |
| **parallel**               | `mclapply()`, `parLapply()` | `mclapply(1:10, f, mc.cores=4)`                            |
| **foreach** + `doParallel` | `%dopar%`                   | `foreach(i=1:10) %dopar% f(i)`                             |
| **future** + `furrr`       | `future_map()`              | `future_map(1:10, f, .options = furrr_options(workers=4))` |
| **BiocParallel**           | `bplapply()`                | `bplapply(X, FUN, BPPARAM=MulticoreParam(4))`              |

Most of these require you to explicitly **enable** parallelization.
Many bioinformatics tools (e.g., **DESeq2**, **MAST**) support this through `BPPARAM = MulticoreParam(n)`.

---

## 4. Quick Test: Is Multicore Working?

```
library(parallel)
system.time(mclapply(1:8, function(i) Sys.sleep(1), mc.cores = 8))
```

- If it finishes in ~1 second → multicore works.
- If it takes ~8 seconds → running serially.

---

## 5. Best Practice on Tufts HPC

When submitting an R job, request the number of cores you plan to use:

```
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
```

Then in your R script:

```
library(parallel)
numCores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
mclapply(1:10, function(x) my_function(x), mc.cores = numCores)
```

This ensures R uses only the cores allocated by SLURM — keeping your job efficient and fair to other users.

---

## Summary

| Step               | Action                              |
| ------------------ | ----------------------------------- |
| Estimate needs     | Start with 4–8 cores                |
| Detect cores       | `parallel::detectCores()`           |
| Use multithreading | Via BLAS or `mclapply()` etc.       |
| Control cores      | `Sys.getenv("SLURM_CPUS_PER_TASK")` |
| Optimize usage     | Check runtime and CPU utilization   |

---
