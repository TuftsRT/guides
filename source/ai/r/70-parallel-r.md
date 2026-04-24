---
tags: r tutorial parallel-computing performance
---

# Parallel Processing in R

**Author:** Kyle Monahan

---

**Learning objectives:** By the end of this tutorial you will understand why and when to parallelize R code, know the difference between FORK and PSOCK cluster types, use `foreach` with `%do%` (sequential) and `%dopar%` (parallel), set up a PSOCK cluster with `doParallel`, benchmark and profile R code, and recognize when parallelization helps, and when it does not.

---

## Introduction

You may have had the experience of waiting too long for a process to run during your R session. One common place this happens is in `for` loops. In this tutorial, we discuss how to speed up loops and other computations in R using parallel processing.

**Packages used:**

- `parallel`: base R parallel library (good for general use; requires manual environment export)
- `doParallel`: a `foreach` backend that handles environment export automatically
- `foreach`: required for `foreach` loops
- `parallelly`: helper functions for parallel work

```r
library(parallel)
library(doParallel)
library(foreach)
library(parallelly)
```

---

## What is a `for` Loop?

Suppose we want to compute the square root of the numbers 1 through 10 and store the results.

```r
x <- vector()           # Create an empty vector outside the loop
for (i in 1:10) {
  x[i] <- sqrt(i)       # Compute square root and store in x
}
x
```

This works fine for small tasks, but by default R uses only a **single core**. Check how many cores your machine has:

```r
parallel::detectCores()
```

If you have 8 or more cores but are using only one, you may be leaving a lot of performance on the table for computationally intensive tasks.

---

## Introducing `foreach`

The `foreach` package provides an alternative to `for` loops that can easily be parallelized. The `%do%` operator runs the loop sequentially on a single core:

```r
x <- foreach::foreach(i = 1:10) %do% {
  sqrt(i)
}
x
```

The result is a list. To combine the outputs into a single vector, use the `.combine` argument:

```r
x <- foreach(i = 1:10, .combine = 'c') %do% {
  sqrt(i)
}
x
```

---

## Running in Parallel with `%dopar%`

To run the loop in parallel, switch `%do%` to `%dopar%`:

```r
x <- foreach(i = 1:10, .combine = 'c') %dopar% {
  sqrt(i)
}
```

This will produce an error! Before using `%dopar%`, you must register a parallel backend (that is, tell R how many cores to use and how they should communicate).

### FORK vs. PSOCK Clusters

There are two main cluster types:

**FORK:**

- Only available on Unix/Linux and macOS
- Workers share memory with the main process (no need to copy the environment)
- Generally faster than PSOCK

**PSOCK (Parallel Socket Cluster):**

- Available on all platforms, including Windows
- Each worker gets its own copy of the environment (can reduce efficiency by ~50% for memory-heavy tasks)
- Recommended for portability and HPC environments

We will use **PSOCK** so the code works on any platform.

### Setting Up a PSOCK Cluster

```r
n.cores    <- parallelly::availableCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
print(my.cluster)
```

Register the cluster as the `doParallel` backend:

```r
doParallel::registerDoParallel(cl = my.cluster)

foreach::getDoParRegistered()  # Should return TRUE
foreach::getDoParWorkers()     # Reports the number of registered workers
```

---

## Comparing Sequential vs. Parallel Execution

Run the loop in parallel:

```r
t1 <- proc.time()

x <- foreach(i = 1:10, .combine = 'c') %dopar% {
  sqrt(i)
}
x

proc.time() - t1
```

And sequentially:

```r
t1 <- proc.time()

x <- foreach(i = 1:10, .combine = 'c') %do% {
  sqrt(i)
}
x

proc.time() - t1
```

The times are roughly the same, or the parallel version may even be *slower*. Computing 10 square roots is trivially fast. The overhead of setting up the parallel cluster and distributing work across cores costs more time than the computation itself.

> **Key insight:** Parallelization only helps when the computation inside each iteration is substantial. For cheap operations, the overhead dominates.

---

## Profiling R Code

Use `Rprof()` to profile a code block:

```r
Rprof(filename = "Profile1.out", line.profiling = TRUE, memory.profiling = TRUE)

x <- foreach(i = 1:10, .combine = 'c') %dopar% { sqrt(i) }
print(x)

Rprof(NULL)
summaryRprof("Profile1.out")
```

You can also use the `profr` package for more detailed profiling output:

```r
library(profr)
p <- profr(
  foreach(i = 1:10, .combine = 'c') %dopar% { sqrt(i) }
)
print(p)
```

RStudio also has a built-in profiler: **Profile > Start Profiling**.

---

## Benchmarking Multiple Approaches

The `rbenchmark` package makes it easy to compare the timing of multiple approaches across many replications:

```r
library(rbenchmark)
library(ggplot2)

set.seed(42)
n.cores <- parallelly::availableCores() - 1
cl      <- makeCluster(n.cores)
registerDoParallel(cl)

FUN <- function(x) { round(sqrt(x), 4) }
a   <- lapply(1:10, function(i) i)

test1 <- benchmark(
  "lapply"        = lapply(1:10, FUN = FUN),
  "For loop"      = for (i in 1:10) { FUN(i) },
  "Foreach dopar" = foreach(i = 1:10) %dopar% FUN(i),
  "Foreach do"    = foreach(i = 1:10) %do% FUN(i),
  "parLapply"     = parLapply(cl = cl, X = a, fun = FUN),
  "parSapply"     = parSapply(cl = cl, X = a, FUN = FUN),
  columns         = c('test', 'elapsed', 'replications'),
  replications    = c(100, 200, 500, 1000)
)

ggplot() +
  geom_line(aes(x = replications, y = elapsed, colour = test), data = test1)
```

For this simple task, `for` loops and `lapply` perform best. The parallel overhead makes `%dopar%` slower.

---

## When Parallelization Does Help: Model Fitting

Parallelization pays off for computationally intensive tasks. Compare approaches for bootstrap resampling of a logistic regression:

```r
FUN <- function(i) {
  ind     <- sample(100, 100, replace = TRUE)
  result1 <- glm(Species ~ Sepal.Length,
                 family = binomial(logit),
                 data   = iris[ind, ])
  coefficients(result1)
}

test5 <- benchmark(
  "lapply"        = lapply(1:10, FUN = FUN),
  "For loop"      = for (i in 1:10) { FUN(i) },
  "Foreach dopar" = foreach(i = 1:10) %dopar% FUN(i),
  "Foreach do"    = foreach(i = 1:10) %do% FUN(i),
  "parLapply"     = parLapply(cl = cl, X = a, fun = FUN),
  "parSapply"     = parSapply(cl = cl, X = a, FUN = FUN),
  columns         = c('test', 'elapsed', 'replications'),
  replications    = c(100, 200, 500, 100)
)

ggplot() +
  geom_line(aes(x = replications, y = elapsed, colour = test), data = test5)
```

Now `parLapply` and `parSapply` are the winners. The computation per iteration is expensive enough to justify the parallel overhead.

---

## Stopping the Cluster

Always stop the cluster when finished to release cores for other processes:

```r
parallel::stopCluster(cl = my.cluster)
```

To release cores without stopping, register a sequential backend:

```r
registerDoSEQ()
```

---

## Parallelizing Mixed Models with `lme4`

Some packages like `lme4` can themselves leverage multiple cores. Compare single-core vs. multi-core model fitting:

**Single core:**

```r
set.seed(42)
ptm        <- proc.time()
cake_model <- lme4::lmer(angle ~ recipe * temperature + (1 | recipe:replicate),
                         data = lme4::cake, REML = FALSE)
proc.time() - ptm
```

**Multi-core (register a parallel backend first):**

```r
doParallel::registerDoParallel(parallel::detectCores() - 1)

ptm   <- proc.time()
model <- lme4::lmer(formula = cake_model,
                    data    = lme4::cake,
                    control = lme4::lmerControl(optimizer = "nloptwrap"))
proc.time() - ptm
```

Use `lme4::allFit()` to check for the optimal model fit across multiple optimizers in parallel:

```r
require(optimx)
require(dfoptim)
nCPU      <- detectCores() - 1
ptm       <- proc.time()
cake_fit  <- lme4::allFit(object  = model,
                          data    = lme4::cake,
                          verbose = TRUE,
                          parallel = 'multicore',
                          ncpus    = nCPU)
proc.time() - ptm
```

```r
cake_fit$bobyqa
```

Release the cores:

```r
registerDoSEQ()
```

---

## Summary

| Concept                     | Key takeaway                                                             |
| --------------------------- | ------------------------------------------------------------------------ |
| **When to parallelize**     | Only when each iteration is computationally expensive                    |
| **FORK vs. PSOCK**          | FORK is faster on Unix/Mac; PSOCK is portable and works everywhere       |
| **`%do%` vs. `%dopar%`**    | `%do%` is sequential; `%dopar%` distributes work across registered cores |
| **`parLapply`/`parSapply`** | Often the fastest for parallel apply-style operations                    |
| **Always stop the cluster** | Call `stopCluster()` when finished                                       |
| **Benchmarking**            | Use `rbenchmark` or `Rprof()` to measure what actually speeds up         |

For support with parallel computing at Tufts, email **datalab-support@elist.tufts.edu**.
