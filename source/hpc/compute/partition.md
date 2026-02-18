---
title: "Partitions and Limits"
---

# Tufts HPC Partitions and Limits


Tufts HPC cluster resources are grouped into different partitions. A partition is a logical collections of nodes that comprise different hardware resources and limits based on **functionality** and **priority** levels.

```{warning}
Research Lab members may have access to additional resources and different limits. See "Contribute Nodes" below.
```

### General limits

**This limit is subject to change to best utilize cluster resources.**

Current Cluster Resource Limits:

> - **Public Partitions** (batch+largemem+gpu)
>
> - CPU: 250 cores
>
>   RAM: 5000 GB
>
>   GPU: 10
>
> - **Preempt Partition** (preempt)
>
> - CPU: 2000 cores
>
>   RAM: 8000 GB
>
>   GPU: 20





## Partitions

**Public Partitions:**

All users have equal access to the following public partitions. Job priorities are under the governance of Slurm Fairshare algorithm.

- **batch**\*: The default partition for standard jobs that do not require any special hardware or configurations. CPU only. Provides memory (RAM) up to 500GB.
- **gpu**: Designated for jobs that require GPU resources. No CPU only jobs allowed.
- **interactive**: Intended for CPU only interactive jobs, which are typically shorter in duration and smaller in resource requirement.
- **largemem**: For CPU only jobs which require up to 1TB of memory (RAM).
- **preempt**: Contains most resources on HPC cluster (CPU and GPU, public and contrib nodes). Jobs submitted to preempt partition has lower priority and can be preempted by contrib node owners' higher priority jobs.

To get a full inventory of available resources and node specs, go to [**OnDemand**](https://ondemand-prod.pax.tufts.edu) `Cluster` --> `System Status`

From command line, use the following command to check what partitions you have access to:

```
$ sinfo
```

## Quality of Service (QOS)



## Time limit

Each partition has a specific time limit.

```
PARTITION       TIMELIMIT
batch*          2-00:00:00
gpu             7-00:00:00
interactive     4:00:00
largemem        7-00:00:00
preempt         7-00:00:00
```

- **preempt** - Be aware, `preempt` partition consists of most of the nodes on the cluster, including **contrib nodes** from different research labs. When submitting jobs to preempt partition, you acknowledge that your jobs are taking the **risk** of being preempted by higher priority jobs. In that case, you will simply have to resubmit your jobs.

# Lab Partitions

Some research labs have dedicated nodes available in the HPC Cluster through our contrib node program. These are accessed using a partition name for each lab.  You can see this name by running the `sinfo` command.

We always reccomend also selecting a public parition in case your lab resources are fully utilized.  Multiple partitions can be specified as a comma seperated list.

`sbatch -p labpartition,public`

```{warning}
Lab partitions may have different resource limits that are more or less restrictive than the defaults above.