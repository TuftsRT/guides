---
tags: python parallel-computing hpc data-science machine-learning
---

# Parallel Computing in Python for Data Science

**Author:** Research Technology, TTS

---

## Learning Objectives

By the end of this tutorial, you will be able to:

1. Explain why parallel processing is important for data science workflows and identify which type of parallelism fits a given task.
1. Set up a reproducible Python environment for parallel computing on the Tufts HPC cluster.
1. Use `multiprocessing`, `joblib`, `concurrent.futures`, and `Dask` to parallelize machine learning and NLP tasks.
1. Apply chunked reading, memory-mapped arrays, and distributed array computation to work with datasets that exceed available RAM.
1. Submit parallel Python jobs to SLURM, including array jobs and jobs scheduled with `--begin`.
1. Offload computation to a GPU using `CuPy` for array operations and `PyTorch` for deep learning and NLP model training.

---

## Overview

Modern data science workflows involving machine learning and natural language processing frequently require training models across large corpora, searching hyperparameter spaces, and running inference over millions of records. Python provides several libraries for parallel computing, from single-machine multiprocessing to distributed cluster computing, making it well suited for HPC environments like the Tufts HPC cluster.

This tutorial covers:

- **Why parallel processing is essential in data science**.
- **Key Python packages** for parallel data science workflows.
- **Code examples to set up and use parallel processing**.
- **Working with datasets larger than memory**.
- **Considerations and best practices**.
- **Applications of parallel processing in machine learning and NLP**.

---

## Why Parallel Processing is Important in Data Science

Data science applications often involve tasks such as:

- Training and evaluating models across large hyperparameter grids or cross-validation folds.
- Tokenizing, embedding, and featurizing large text corpora for NLP pipelines.
- Running inference at scale across millions of records or documents.

These tasks can be resource-intensive, requiring significant memory and processing power. Parallel processing distributes work across multiple cores or nodes, reducing overall computation time.

---

## Key Python Packages for Parallel Processing in Data Science

Several Python packages facilitate parallel processing, each suited for different types of workflows and computing environments:

1. **`multiprocessing`**: Python's built-in module for process-based parallelism on a single machine.
1. **`joblib`**: A high-level library for easy parallelization of loops and pipelines, widely used in scikit-learn and data science tooling.
1. **`concurrent.futures`**: A clean, modern interface for thread- and process-based parallelism in the standard library.
1. **`Dask`**: A flexible library for parallel and distributed computing that scales from a laptop to a SLURM cluster, with support for out-of-memory datasets.

---

## Setting Up Parallel Processing in Python

### Environment Setup

We assume you know how to access the Tufts HPC cluster. If you don't, please see these resources:

- [Accessing and Using the HPC Cluster](https://rtguides.it.tufts.edu/hpc/access/index.html)
- [HPC Examples and Tutorials](https://rtguides.it.tufts.edu/hpc/examples/index.html)

Before running parallel jobs, set up a reproducible Python environment. The recommended approach on the Tufts HPC cluster is to use `miniforge`, then redirect conda environments and packages to your lab's research storage to avoid filling your home directory quota. See the [Conda Environments guide](https://rtguides.it.tufts.edu/hpc/application/10-condaenv.html) for full instructions.

**Step 1: Start an interactive session on a compute node.** Never install packages on the login node.

```bash
srun --pty -t 2:00:00 -p batch --mem=8g -N 1 -n 4 bash
```

**Step 2: Load the recommended module.**

```bash
module load miniforge/25.3.0
```

**Step 3: Configure conda to store environments and packages in your lab storage** (only needed once). Replace `XXXXlab` with your actual lab group name:

```bash
conda config --append envs_dirs /cluster/tufts/XXXXlab/$USER/condaenv/
conda config --append pkgs_dirs /cluster/tufts/XXXXlab/$USER/condapkg/
conda config --add channels conda-forge
```

**Step 4: Create and activate your environment.**

```bash
conda create -n ds_parallel python=3.11 numpy pandas scikit-learn dask joblib -y
conda activate ds_parallel
```

> **Note:** If a package is available via conda (search [Anaconda Cloud](https://anaconda.org/search)), prefer `conda install` over `pip install`. Mixing conda and pip can corrupt environments. Only use `pip` for packages not available through conda channels.

There are many packages to implement parallel processing in Python. If you are just getting started, we suggest you start with `Dask` (section 4, below).

### 1. Using `multiprocessing`

The `multiprocessing` module spawns independent worker processes, bypassing Python's Global Interpreter Lock (GIL) and making it suitable for CPU-bound tasks.

```python
from multiprocessing import Pool
import os


def train_on_fold(fold_idx):
    """Train and evaluate a model on one CV fold."""
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import roc_auc_score
    import numpy as np

    rng = np.random.default_rng(fold_idx)
    # In practice, load your actual fold data here
    X_train = rng.standard_normal((1000, 50))
    y_train = rng.integers(0, 2, 1000)
    X_val = rng.standard_normal((200, 50))
    y_val = rng.integers(0, 2, 200)

    clf = RandomForestClassifier(n_estimators=200, random_state=fold_idx)
    clf.fit(X_train, y_train)
    return roc_auc_score(y_val, clf.predict_proba(X_val)[:, 1])


if __name__ == "__main__":
    n_cores = os.cpu_count() - 1
    with Pool(processes=n_cores) as pool:
        fold_aucs = pool.map(train_on_fold, range(10))
    print(f"Mean AUC: {sum(fold_aucs) / len(fold_aucs):.4f}")
```

The `if __name__ == "__main__":` guard is required when using `multiprocessing` on Linux to prevent worker processes from re-importing and re-running the script.

### 2. Using `joblib`

`joblib` provides a simple `Parallel` + `delayed` interface that is especially convenient for parallelizing loops with minimal boilerplate. It is the backend used by scikit-learn and is ideal for hyperparameter search and feature engineering pipelines.

```python
from joblib import Parallel, delayed
import numpy as np
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score

param_grid = [
    {"n_estimators": n, "max_depth": d, "learning_rate": lr}
    for n in [100, 300, 500]
    for d in [3, 5]
    for lr in [0.01, 0.1]
]


def evaluate_config(params, X_train, y_train, X_val, y_val):
    clf = GradientBoostingClassifier(**params, random_state=42)
    clf.fit(X_train, y_train)
    return {
        "params": params,
        "auc": roc_auc_score(y_val, clf.predict_proba(X_val)[:, 1]),
    }


# n_jobs=-2 leaves one core free
results = Parallel(n_jobs=-2, verbose=5)(
    delayed(evaluate_config)(p, X_train, y_train, X_val, y_val) for p in param_grid
)

best = max(results, key=lambda r: r["auc"])
print(f"Best AUC: {best['auc']:.4f} | Params: {best['params']}")
```

`verbose=5` prints progress updates, which is useful for long-running searches.

### 3. Using `concurrent.futures`

`concurrent.futures` offers a clean, future-based API. Use `ProcessPoolExecutor` for CPU-bound work and `ThreadPoolExecutor` for I/O-bound work such as downloading datasets or reading many files in parallel.

```python
from concurrent.futures import ProcessPoolExecutor, as_completed
import os


def preprocess_document(filepath):
    """Tokenize and clean a single document file."""
    with open(filepath) as f:
        text = f.read().lower()
    tokens = [w for w in text.split() if w.isalpha()]
    return filepath, tokens


doc_paths = [f"corpus/doc_{i}.txt" for i in range(1000)]

with ProcessPoolExecutor(max_workers=os.cpu_count() - 1) as executor:
    futures = {executor.submit(preprocess_document, p): p for p in doc_paths}

    token_store = {}
    for future in as_completed(futures):
        path = futures[future]
        try:
            filepath, tokens = future.result()
            token_store[filepath] = tokens
        except Exception as e:
            print(f"Error processing {path}: {e}")
```

### 4. Using `Dask` for Scalable Parallelism

`Dask` provides parallelism that scales from a single machine to a SLURM cluster without changing your code. It is particularly well suited for dataframe and array operations on data too large to fit in memory, which is common in large-scale NLP feature pipelines and model evaluation workflows.

#### Single-Machine Dask

```python
import dask.dataframe as dd
from dask.distributed import Client

# Start a local cluster (auto-detects available cores)
client = Client()
print(client.dashboard_link)  # open in browser to monitor progress

# Read a large feature CSV without loading it into RAM
df = dd.read_csv("features_large.csv")

# Operations are lazy -- nothing is computed until .compute() is called
high_conf = df[df["confidence"] > 0.9][["doc_id", "label", "confidence"]]
result = high_conf.compute()
print(result.head())

client.close()
```

#### Scaling Dask to a SLURM Cluster

Using `dask-jobqueue`, you can launch Dask workers as SLURM jobs directly from Python. Install it first:

```bash
conda install -c conda-forge dask-jobqueue
```

```python
from dask_jobqueue import SLURMCluster
from dask.distributed import Client

cluster = SLURMCluster(
    queue="batch,preempt",  # use preempt for faster queue times when possible
    cores=8,  # cores per worker job
    memory="32GB",  # memory per worker job
    walltime="04:00:00",
    job_extra_directives=["--output=dask_worker_%j.log"],
)

# Scale to 4 worker jobs (4 x 8 = 32 total cores)
cluster.scale(jobs=4)
client = Client(cluster)

# Run a distributed TF-IDF computation across the cluster
import dask.array as da
import dask.bag as db

corpus = db.read_text("corpus/*.txt")
word_counts = corpus.map(str.split).flatten().frequencies().compute()

client.close()
cluster.close()
```

---

## Working with Datasets Larger Than Memory in Python

Large data science datasets (web-scraped corpora, high-dimensional embedding stores, large feature matrices) frequently exceed available RAM. The strategies below allow you to work with such datasets without loading them fully into memory.

### Understanding the Three Strategies

| Strategy                          | How it works                                            | Best for                                                    |
| --------------------------------- | ------------------------------------------------------- | ----------------------------------------------------------- |
| **Chunked / streaming reads**     | Data is read one piece at a time                        | Large tabular files (feature CSVs, log data, text datasets) |
| **Memory-mapped arrays**          | OS pages in data on demand from disk                    | Large numeric arrays needing random access                  |
| **Distributed array computation** | Lazy operations evaluated block-by-block across workers | Any dataset; integrates with Dask parallelism               |

### Strategy 1: Chunked Reads with `pandas`

```python
import pandas as pd

chunk_size = 500_000
results = []

for chunk in pd.read_csv("large_features.csv", chunksize=chunk_size):
    # Filter and aggregate each chunk before accumulating
    sig = chunk[chunk["confidence"] > 0.9]
    results.append(sig)

filtered = pd.concat(results, ignore_index=True)
print(filtered.shape)
```

### Strategy 2: On-Disk Arrays with `h5py`

HDF5 is a common format for storing pre-computed embeddings and large numeric arrays. `h5py` exposes HDF5 datasets as NumPy-like arrays that are read from disk only when sliced.

```python
import h5py
import numpy as np

with h5py.File("document_embeddings.h5", "r") as f:
    f.visit(print)  # inspect structure

    embeddings = f["embeddings"]
    print(f"Embedding matrix shape: {embeddings.shape}")

    # Read only the first 10,000 document embeddings into RAM
    subset = embeddings[:10_000, :]
```

For datasets stored in the Hugging Face `datasets` format, use memory-mapped Arrow files:

```python
from datasets import load_from_disk

# Loaded in memory-mapped mode by default -- only accessed data is read
ds = load_from_disk("large_text_dataset")
sample = ds.select(range(1000))
```

### Strategy 3: Distributed Arrays with `Dask`

Dask arrays apply lazy operations block-by-block across available workers, making them ideal for large embedding matrices or feature stores.

```python
import dask.array as da
import h5py

# Open an HDF5 embedding matrix as a Dask array (no data loaded yet)
f = h5py.File("document_embeddings.h5", "r")
embeddings = da.from_array(f["embeddings"], chunks=(10_000, 768))

# Compute per-dimension means block-by-block
dim_means = embeddings.mean(axis=0).compute()
print(f"Embedding dimension means shape: {dim_means.shape}")

f.close()
```

Combine with the SLURM-backed Dask cluster above to distribute computation across multiple nodes.

#### Practical Tips

- **Prefer HDF5 or Parquet over plain CSV** for large datasets. Both support partial reads; Parquet also supports column pruning and predicate pushdown via `pyarrow`.
- **Tune chunk sizes** to fit comfortably in per-worker memory. A useful starting point is `chunk_size = available_memory_per_worker / n_workers / 4`.
- **Avoid implicit `.values` or `np.array()` calls** on Dask arrays; these trigger full materialisation into RAM.
- **Use Hugging Face datasets in memory-mapped mode** by default; avoid `.to_pandas()` on large splits unless you have sufficient RAM.

---

## Considerations and Best Practices

1. **Memory usage**: Each worker process gets its own memory space. Ensure your SLURM allocation accounts for the number of workers times the per-worker footprint, not just the dataset size.
1. **Overhead**: Spawning processes has startup cost. For very short tasks, `joblib` with `prefer="threads"` or simple vectorized NumPy/pandas operations will outperform process-based parallelism.
1. **Error handling**: Errors in worker processes can be silent. Use `concurrent.futures`'s `future.result()` inside a `try/except`, or `joblib`'s `return_as="generator"` to catch per-task failures.
1. **The GIL**: Python's Global Interpreter Lock prevents true thread-based CPU parallelism. Always use process-based backends (`multiprocessing`, `ProcessPoolExecutor`, `joblib` with `backend="loky"`) for CPU-bound work.
1. **Reproducibility**: Set random seeds inside each worker explicitly; workers do not inherit the parent process's random state.
1. **System load**: On shared nodes, leave at least one core free (`n_jobs=-2` in `joblib`, `os.cpu_count() - 1` elsewhere). On SLURM, request only the cores you need.

---

## Applications in Data Science

### Hyperparameter Search

Parallelise a random search over model configurations with `joblib`:

```python
from joblib import Parallel, delayed
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
import numpy as np

rng = np.random.default_rng(42)
param_samples = [
    {
        "n_estimators": rng.integers(100, 1000),
        "max_depth": rng.integers(3, 20),
        "max_features": rng.choice(["sqrt", "log2"]),
    }
    for _ in range(50)
]


def fit_and_score(params):
    clf = RandomForestClassifier(**params, random_state=0, n_jobs=1)
    clf.fit(X_train, y_train)
    return {
        "params": params,
        "auc": roc_auc_score(y_val, clf.predict_proba(X_val)[:, 1]),
    }


results = Parallel(n_jobs=-2, verbose=5)(
    delayed(fit_and_score)(p) for p in param_samples
)
best = max(results, key=lambda r: r["auc"])
print(best)
```

### Cross-Validation

Scikit-learn's `cross_val_score` uses `joblib` under the hood; set `n_jobs` directly:

```python
from sklearn.pipeline import Pipeline
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score

pipeline = Pipeline(
    [
        ("tfidf", TfidfVectorizer(max_features=50_000, ngram_range=(1, 2))),
        ("clf", LogisticRegression(max_iter=1000)),
    ]
)

# n_jobs=-1 distributes folds across all available cores
scores = cross_val_score(pipeline, texts, labels, cv=10, scoring="roc_auc", n_jobs=-1)
print(f"CV AUC: {scores.mean():.4f} (+/- {scores.std():.4f})")
```

### Text Preprocessing and Tokenization

Distribute tokenization of a large corpus across cores with `multiprocessing`:

```python
from multiprocessing import Pool
import os
from pathlib import Path


def tokenize_file(path):
    text = Path(path).read_text()
    return path, [w.lower() for w in text.split() if w.isalpha()]


doc_paths = list(Path("corpus").glob("*.txt"))

with Pool(processes=os.cpu_count() - 1) as pool:
    token_map = dict(pool.map(tokenize_file, doc_paths))

print(f"Tokenized {len(token_map)} documents")
```

### Distributed Inference with a Transformer Model

Run inference over a large document collection by batching and parallelising across cores. For GPU inference, see the GPU section below.

```python
from joblib import Parallel, delayed
from transformers import pipeline


def embed_batch(batch):
    """Run sentence embedding on a batch of texts."""
    embedder = pipeline(
        "feature-extraction", model="sentence-transformers/all-MiniLM-L6-v2"
    )
    return [out[0][0] for out in embedder(batch)]


# Split corpus into batches of 64
batch_size = 64
batches = [texts[i : i + batch_size] for i in range(0, len(texts), batch_size)]

embeddings = Parallel(n_jobs=-2, verbose=10)(delayed(embed_batch)(b) for b in batches)
```

### NLP Model Fine-Tuning

Use `Dask` to preprocess a large text dataset before fine-tuning, keeping memory usage bounded:

```python
import dask.bag as db
from transformers import AutoTokenizer

tokenizer = AutoTokenizer.from_pretrained("bert-base-uncased")

corpus = db.read_text("corpus/*.txt")

tokenized = corpus.map(
    lambda text: tokenizer(text, truncation=True, max_length=512, return_tensors="pt")
).compute()
```

### Bootstrap Confidence Intervals

Parallelise bootstrapping for model evaluation metrics:

```python
from joblib import Parallel, delayed
import numpy as np
from sklearn.metrics import roc_auc_score


def bootstrap_auc(i, y_true, y_prob):
    rng = np.random.default_rng(i)
    idx = rng.integers(0, len(y_true), len(y_true))
    return roc_auc_score(y_true[idx], y_prob[idx])


boot_aucs = Parallel(n_jobs=-2)(
    delayed(bootstrap_auc)(i, y_test, y_prob) for i in range(2000)
)

lo, hi = np.percentile(boot_aucs, [2.5, 97.5])
print(f"AUC 95% CI: {lo:.4f} - {hi:.4f}")
```

---

## Setting Up Jobs on Open OnDemand or SLURM

### Understanding Node Hardware

Before submitting a job, inspect available resources with `hpctools`:

```bash
module load hpctools
hpctools
```

Select `1. Checking Free Resources On Each Node in Given Partition(s)` then `batch` to see a table like:

```
NODELIST     STATE  PARTITION  MEMORY   ALLOCMEM  CPUS(A/I/O/T)
p1cmp018     mix    batch*     248000   244445    16/56/0/72
p1cmp019     mix    batch*     248000   244445    16/56/0/72
d1cmp028     mix    batch*     510000   506880    88/40/0/128
```

`CPUS(A/I/O/T)` shows allocated, idle, other, and total cores.

### SLURM Job Script for a Parallel Python Job

```bash
#!/bin/bash -l
#SBATCH -J ds_parallel                 # job name
#SBATCH --time=00-12:00:00             # requested time (DD-HH:MM:SS)
#SBATCH -p batch,preempt               # run on batch or preempt, whichever is free first
#SBATCH -N 1                           # 1 node (leave at 1 for shared-memory jobs)
#SBATCH -n 8                           # 8 tasks = 8 CPU cores for this job
#SBATCH --mem=64g                      # total RAM for the job
#SBATCH --output=MyJob.%j.%N.out       # standard output, %j=JOBID %N=NodeName
#SBATCH --error=MyJob.%j.%N.err        # standard error
#SBATCH --mail-type=ALL                # email on start, end, and fail
#SBATCH --mail-user=Your_Tufts_Email@tufts.edu

module purge
module load miniforge/25.3.0
source activate ds_parallel

python my_parallel_script.py

conda deactivate
```

Key parameters:

- `-n 8`: the number of CPU cores allocated. Match this to `n_jobs` or `max_workers` in your script.
- `--mem=64g`: total memory for the job. For parallel jobs, budget per-worker memory times the number of workers.
- `--time=00-12:00:00`: wall-clock limit in DD-HH:MM:SS format. Set conservatively for first runs.
- `-p batch,preempt`: listing both partitions lets SLURM start the job on whichever has resources available first.

Submit with:

```bash
sbatch my_parallel_script.sh
```

### SLURM Array Jobs

For independent model runs (one job per fold, dataset, or random seed), SLURM array jobs are more efficient than spawning many Python processes from a single script:

```bash
#!/bin/bash -l
#SBATCH -J cv_array
#SBATCH --time=00-04:00:00
#SBATCH -p batch,preempt
#SBATCH --array=1-10                   # one task per CV fold
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=32g
#SBATCH --output=fold%a.%j.%N.out
#SBATCH --error=fold%a.%j.%N.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Your_Tufts_Email@tufts.edu

module purge
module load miniforge/25.3.0
source activate ds_parallel

python train_fold.py --fold ${SLURM_ARRAY_TASK_ID}

conda deactivate
```

`${SLURM_ARRAY_TASK_ID}` is automatically set to the task index, allowing each job to process a different fold or configuration independently.

---

## Automatic Core Detection

Python's `os.cpu_count()` returns the total number of logical cores on the machine. On SLURM, however, this reflects the node's total core count rather than your allocation. To respect your SLURM allocation:

```python
import os

# Prefer the SLURM-allocated core count; fall back to os.cpu_count()
n_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count()))
print(f"Using {n_workers} workers")
```

Always pass `n_workers` explicitly to `Pool`, `Parallel`, or `ProcessPoolExecutor` rather than relying on auto-detection when running under SLURM.

---

## GPU Computing in Python on HPC

Many data science workloads (training large language models, running transformer inference, computing embeddings at scale) run significantly faster on a GPU. The Tufts HPC cluster provides GPU nodes accessible via SLURM, including A100, L40, L40S, and H200 GPUs. This section covers two complementary tools: `CuPy` for drop-in GPU acceleration of NumPy-style array operations, and `PyTorch` for deep learning and NLP model training.

### Requesting a GPU Node via SLURM

All GPU jobs require the `gpu` partition and a `--gres` directive specifying the GPU type and count. Use `module load cuda/12.9.0` to provide the necessary GPU libraries.

```bash
#!/bin/bash -l
#SBATCH -J gpu_ds
#SBATCH --time=00-08:00:00
#SBATCH -p gpu,preempt                      # gpu or preempt partition
#SBATCH -N 1
#SBATCH -n 4                                # CPU cores to support the GPU workers
#SBATCH --mem=64g
#SBATCH --gres=gpu:a100:1                   # request 1 A100 GPU
#SBATCH --constraint="a100-80G"             # request the 80GB variant
#SBATCH --begin=2026-04-10T09:00:00         # start no earlier than this time (YYYY-MM-DDTHH:MM:SS)
#SBATCH --output=MyJob.%j.%N.out
#SBATCH --error=MyJob.%j.%N.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Your_Tufts_Email@tufts.edu

module purge
module load miniforge/25.3.0
module load cuda/12.9.0
source activate ds_parallel

python my_gpu_script.py

conda deactivate
```

Use `hpctools` to check which GPU nodes are available and what GPU models they carry before submitting.

The `--begin` directive tells SLURM not to start the job before a specified time. It accepts an absolute timestamp (`YYYY-MM-DDTHH:MM:SS` in local cluster time) or relative values such as `now+2hours` or `16:00:00` (today at 4 PM). Until the begin time is reached, the job sits in a pending (`PD`) state with reason `BeginTime`. After that point it enters the normal queue and competes for available resources, so it will not necessarily start exactly at the specified time if the requested node is busy.

---

### CuPy: NumPy on the GPU

`CuPy` is a GPU array library with an interface nearly identical to NumPy. For code that already uses NumPy, switching to CuPy often requires changing only the import and array creation; the rest of the code stays the same.

#### Installation

Install CuPy matching the cluster's CUDA version. Since `cuda/12.9.0` is the current module, use:

```bash
pip install cupy-cuda12x
```

#### Bioinformatics Example: Large Embedding Matrix Operations

GPU acceleration is particularly effective for the dense linear algebra common in NLP, including similarity computation, dimensionality reduction, and clustering on embedding matrices:

```python
import cupy as cp
import numpy as np

# Load a pre-computed document embedding matrix on CPU (docs x dims)
embeddings_cpu = np.load("document_embeddings.npy").astype(np.float32)

# Transfer to GPU
embeddings_gpu = cp.asarray(embeddings_cpu)

# L2-normalise embeddings entirely on GPU
norms = cp.linalg.norm(embeddings_gpu, axis=1, keepdims=True)
normalised = embeddings_gpu / norms

# Compute full cosine similarity matrix (for smaller corpora)
similarity_matrix = cp.dot(normalised, normalised.T)

# Transfer result back to CPU
sim_cpu = cp.asnumpy(similarity_matrix)
print(f"Similarity matrix shape: {sim_cpu.shape}")
```

#### GPU-Accelerated PCA on Embeddings

```python
import cupy as cp
from cupy.linalg import svd

# Center the embedding matrix
centered = embeddings_gpu - embeddings_gpu.mean(axis=0)

# Truncated SVD for PCA (top 50 dimensions)
# For very large matrices, use cuml.decomposition.PCA (RAPIDS) instead
U, S, Vt = svd(centered, full_matrices=False)
pca_coords = U[:, :50] * S[:50]

pca_cpu = cp.asnumpy(pca_coords)
print(f"PCA embedding shape: {pca_cpu.shape}")
```

#### Checking GPU Memory

```python
import cupy as cp

mempool = cp.get_default_memory_pool()
print(f"Used GPU memory:  {mempool.used_bytes() / 1e9:.2f} GB")
print(f"Total GPU memory: {cp.cuda.Device().mem_info[1] / 1e9:.2f} GB")

# Free cached memory between large operations
mempool.free_all_blocks()
```

#### Practical Tips for CuPy

- **Minimize CPU/GPU transfers.** Each `cp.asnumpy()` or `cp.asarray()` call copies data over the PCIe bus. Chain GPU operations together and transfer only final results.
- **Use `float32` instead of `float64`.** GPUs are optimized for 32-bit arithmetic; most NLP and ML applications do not need double precision and will see 2-4x speedups with `float32`.
- **Fall back gracefully.** Wrap CuPy imports so your code still runs on CPU-only nodes:

```python
try:
    import cupy as cp

    xp = cp  # GPU available
except ImportError:
    import numpy as cp

    xp = cp  # Fall back to NumPy on CPU
```

---

### PyTorch: Deep Learning and NLP Model Training

`PyTorch` is the dominant deep learning framework for NLP and data science, powering transformer-based models such as BERT, GPT, and their derivatives, as well as general-purpose neural networks for tabular and time-series data.

#### Installation

Install PyTorch with CUDA support matching `cuda/12.9.0`:

```bash
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu121
```

#### Device Management

```python
import torch

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

if torch.cuda.is_available():
    print(f"GPU: {torch.cuda.get_device_name(0)}")
    print(f"Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
```

#### Fine-Tuning a Text Classifier with Hugging Face and PyTorch

The most common NLP GPU workload is fine-tuning a pre-trained transformer for text classification:

```python
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from torch.optim import AdamW

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

model_name = "distilbert-base-uncased"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModelForSequenceClassification.from_pretrained(model_name, num_labels=2).to(
    device
)


# Tokenize your dataset
def tokenize(texts, labels, max_length=128):
    enc = tokenizer(
        texts, truncation=True, padding=True, max_length=max_length, return_tensors="pt"
    )
    enc["labels"] = torch.tensor(labels)
    return enc


# Assume texts_train and labels_train are lists
train_enc = tokenize(texts_train, labels_train)
loader = DataLoader(train_enc, batch_size=32, shuffle=True)

optimizer = AdamW(model.parameters(), lr=2e-5)

for epoch in range(3):
    model.train()
    total_loss = 0.0
    for batch in loader:
        batch = {k: v.to(device) for k, v in batch.items()}
        outputs = model(**batch)
        loss = outputs.loss
        loss.backward()
        optimizer.step()
        optimizer.zero_grad()
        total_loss += loss.item()
    print(f"Epoch {epoch+1} | Loss: {total_loss / len(loader):.4f}")

# Save the fine-tuned model
model.save_pretrained("fine_tuned_classifier")
tokenizer.save_pretrained("fine_tuned_classifier")
```

#### Efficient Inference at Scale

For large-scale inference (e.g., classifying millions of documents), use `model.eval()`, `torch.no_grad()`, and batch your inputs:

```python
from transformers import pipeline

classifier = pipeline(
    "text-classification", model="fine_tuned_classifier", device=0
)  # device=0 uses the first GPU

batch_size = 128
predictions = []
for i in range(0, len(texts), batch_size):
    batch = texts[i : i + batch_size]
    predictions.extend(classifier(batch))
```

#### Multi-GPU Training with `DataParallel`

If your SLURM allocation includes multiple GPUs, PyTorch can distribute batches across them with minimal code changes. First, update your SLURM header:

```bash
#SBATCH --gres=gpu:a100:4   # request 4 A100 GPUs
```

Then wrap your model before moving it to the device:

```python
import torch
import torch.nn as nn

device = torch.device("cuda")
model = AutoModelForSequenceClassification.from_pretrained(model_name, num_labels=2)

if torch.cuda.device_count() > 1:
    print(f"Using {torch.cuda.device_count()} GPUs")
    model = nn.DataParallel(model)

model = model.to(device)
```

#### Practical Tips for PyTorch

- **Always call `model.eval()` and `torch.no_grad()` during inference.** This disables dropout and gradient tracking, reducing memory use and speeding up forward passes.
- **Use mixed precision training** for large models. `torch.cuda.amp.autocast()` uses float16 for most operations, roughly halving memory usage and doubling throughput on modern GPUs:

```python
from torch.cuda.amp import autocast, GradScaler

scaler = GradScaler()
for batch in loader:
    batch = {k: v.to(device) for k, v in batch.items()}
    optimizer.zero_grad()
    with autocast():
        outputs = model(**batch)
        loss = outputs.loss
    scaler.scale(loss).backward()
    scaler.step(optimizer)
    scaler.update()
```

- **Save and resume checkpoints** for long training runs. SLURM jobs can be preempted; saving periodically prevents losing hours of training:

```python
# Save
torch.save(
    {"epoch": epoch, "model": model.state_dict(), "optimizer": optimizer.state_dict()},
    "checkpoint.pt",
)

# Resume
ckpt = torch.load("checkpoint.pt")
model.load_state_dict(ckpt["model"])
optimizer.load_state_dict(ckpt["optimizer"])
start_epoch = ckpt["epoch"] + 1
```

- **Monitor GPU utilization** during a job with `nvidia-smi` in a separate terminal or via the Open OnDemand dashboard. Utilization consistently below 50% suggests the data loading pipeline (not the GPU) is the bottleneck; try increasing `DataLoader` workers: `DataLoader(..., num_workers=4, pin_memory=True)`.
