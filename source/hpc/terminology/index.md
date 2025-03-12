# HPC Terminology

## What is a "Cluster"

- A computer cluster is a set of loosely or tightly **connected computers (Nodes)** that work together so that, in many respects, they can be viewed as a single system. Computer clusters have each computing unit set to perform the same/similar tasks, controlled and scheduled by [**software**](../slurm/index).

```{gallery-grid}
---
grid-columns: 1
---
- header: "{fas}`book;pst-color-primary` CPU vs GPU"
  content: "Difference between CPUs and GPUs"
  link: "cpugpu.html"

- header: "{fas}`book;pst-color-primary` Cores vs Nodes"
  content: "Common terms of cores and nodes used when requesting resources"
  link: "cores-nodes.html"

- header: "{fas}`book;pst-color-primary` Memory vs Storage"
  content: "What is memory or storage? How to request them"
  link: "memory-storage.html"

```

## What to Expect on Tufts HPC Cluster?

- Each of the CPU cores on the HPC cluster could be less powerful than the ones on your laptop. The CPU clock speed on Tufts HPC cluster is ranged between 2.1GHz - 3.1 GHz.
- There are a lot more CPU cores on a single node on the cluster than what's available on your local desktop or laptop.
- If you don't have a framework that can utilize multiple GPUs, DO NOT request more than one GPU.
- The more resources (Memory, CPU cores, GPUs) you request, the longer your job might have to wait in the queue.
- Cluster VAST storage is fast, but may not be as fast as your local SSD.
- When you first login to the HPC, you will be on a **login node**, while computation needs to be completed on **compute node**. You know you are on a login node when you see this via the shell: `your_utln@login-prod-01`. See [here](../access/20-cli) for details on how to access the HPC and using these nodes.
