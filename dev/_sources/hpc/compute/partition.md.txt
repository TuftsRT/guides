---
title: Partitions and Limits
---

# Tufts HPC Partitions and Limits

Tufts HPC cluster resources are grouped into different partitions. A partition is a logical collections of nodes that comprise different hardware resources and limits based on functionality, ownership and priority levels.

```{warning}
Research Lab members may have access to additional resources and different limits. See ["Lab Partitions"](#lab-partitions) section below.
```

Up-to-date system information, real-time status, and resource availability can be viewed with [OnDemand System Status](https://ondemand-p01.pax.tufts.edu/pun/sys/dashboard/system-status) app.

### General limits

In general the HPC resources in the public partitions are available to reseachers on a first come, first serve, basis with jobs submitted first, getting access to the next available resource that matches their request(s). However to fairly balance usage across the community heavy users may see their jobs wait if the cluster is fully utilized, with lower usage users getting priority.

A limit is placed on the total resources a single user can have allocated at any one time. The current limits are shown here.

**Public Partitions** (batch + gpu)

- CPU: 250 cores
- RAM: 5000 GB
- GPU: 10

**Preempt Partition** (preempt)

- CPU: 1000 cores
- RAM: 8000 GB
- GPU: 20

Additionally, each user is allowed to request a single interactive job which will have priority over non interactive jobs.

<small>These limits are subject to change to best optimize the utilization of the cluster resources.</small>

## Partitions

**Public Partitions:**

All users have equal access to the following public partitions. Job priorities are under the governance of Slurm Fairshare algorithm.

- **batch**\*: The default partition for standard jobs that do not require any special hardware or configurations. CPU only. Provides memory (RAM) up to 500GB.
- **gpu**: Designated for jobs that require GPU resources. No CPU only jobs allowed.
- **preempt**: Contains most resources on HPC cluster (CPU and GPU, public and contrib nodes). Jobs submitted to preempt partition has lower priority and can be preempted by contrib node owners' higher priority jobs.

The mpi, largemem, and interactive partitions have been retired. Use the batch or gpu partitions instead.

To get a full inventory of specific available resources and node specs, go to [**OnDemand**](https://ondemand-prod.pax.tufts.edu) `Cluster` --> `System Status`

From command line, use the following command to check what partitions you have access to:

```
$ sinfo
```

## Quality of Service (QOS)

The cluster utilizes Slurm QOS to manage special cases and exceptions to the default resource and time limits.

Common QOSs available are

- Interactive - Single Job for up to 4 hours with 1 GPU,

## Time limit

Each partition has a time limit that a jobs runtime cannot exceed. Most use cases should be able to complete within these time ranges, but we understand their are exception. Please contact us if you need access to longer runtimes.

```
PARTITION       TIMELIMIT
batch*          2-00:00:00
gpu             2-00:00:00
preempt         2-00:00:00
```

- **preempt** - Be aware, `preempt` partition consists of most of the nodes on the cluster, including contrib nodes from different research labs. When submitting jobs to preempt partition, you acknowledge that your jobs are taking the risk of being preempted by higher priority jobs. In that case, you will simply have to resubmit your jobs.

# Lab Partitions

Some research labs have dedicated nodes available in the HPC Cluster through our [contrib node](../policy/contribute-nodes) program. These are accessed using a partition name for each lab. You can see this name by running the `sinfo` command.

We always recommend also selecting a public partition in case your lab resources are fully utilized. Multiple partitions can be specified as a comma separated list.

`sbatch -p labpartition,batch` or `sbatch -p labpartition,gpu`

```{warning}
Lab partitions may have different resource limits that are more or less restrictive than the defaults above.
```
