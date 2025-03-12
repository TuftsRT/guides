# Tufts HPC Partitions

% This is where QOS and Reservations will also go

- Tufts HPC cluster resources are grouped into different partitions. A partition is a logical collections of nodes that comprise different hardware resources and limits based on **functionality** and **priority** levels.

## Partitions

**Public Partitions:**

All users have equal access to the following public partitions. Job priorities are under the governance of Slurm Fairshare algorithm.

- **batch**\*: The default partition for standard jobs that do not require any special hardware or configurations. CPU only. Provides memory (RAM) up to 500GB.
- **gpu**: Designated for jobs that require GPU resources. No CPU only jobs allowed.
- **interactive**: Intended for CPU only interactive jobs, which are typically shorter in duration and smaller in resource requirement.
- **largemem**: For CPU only jobs which require up to 1TB of memory (RAM).
- **mpi**: Duplicate set of CPU resources as batch partition. Mainly for jobs that use the Message Passing Interface (MPI) for parallel computing across multiple nodes.
- **preempt**: Consists most resources on HPC cluster (CPU and GPU, public and contrib nodes). Jobs submitted to preempt partition has lower priority and can be preempted by contrib node owners' higher priority jobs.

To get a full inventory of available resources and node specs, go to [**OnDemand**](https://ondemand.pax.tufts.edu) `Misc` --> `Inventory`

From command line, use the following command to check what partitions you have access to:

```
$ sinfo
```

**Lab Partitions:**

If you have special access to research groups' contrib node, the lab partition will be listed in `sinfo` output in addition to the public partitions mentioned above.

## Time limit

Each partition has a specific time limit.

```
PARTITION       TIMELIMIT
batch*          7-00:00:00
gpu             7-00:00:00
interactive     4:00:00
largemem        7-00:00:00
mpi             7-00:00:00
preempt         7-00:00:00
```

- **preempt** - Be aware, `preempt` partition consists of most of the nodes on the cluster, including **contrib nodes** from different research labs. When submitting jobs to preempt partition, you acknowledge that your jobs are taking the **risk** of being preempted by higher priority jobs. In that case, you will simply have to resubmit your jobs.
