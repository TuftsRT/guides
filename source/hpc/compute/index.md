# Computing Resources

The Tufts HPC Cluster is comprised of various type of CPU and GPU enabled compute nodes. This table provides a
general idea of what is currently available for use.

| Number | Type  | CPU         | Memory | GPU         |
| :----- | :---- | :---------- | :----- | :---------- |
| 40     | CPU   | 2x 32 Cores | 512GB  |             |
| 45     | CPU   | 2x 18 Cores | 256GB+ |             |
| 2      | CPU   | 2x 18 Cores | 1TB    |             |
| 8      | CPU   | 2x 24 Cores | 756GB  |             |
| 8      | CPU   | 2x 32 Cores | 256GB  |             |
| 5      | GPU   | 2x 16 Cores | 756GB  | 8x a100-40G |
| 5      | GPU   | 2x 32 Cores | 1TB    | 8x a100-80G |
| 20     | GPU   | 2x 20 Cores | 92GB   | 4x t4       |
| 9      | GPU   | 2x 32 Cores | 256GB  | 4x l40      |
| 7      | GPU   | 2x 32 Cores | 512GB  | 4x l40s     |
| 4      | GPU   | 2x 48 Cores | 1.5TB  | 8x H200     |
| 25     | MIXED | MIXED       | MIXED  | MIXED       |

Some systems are only available to all users on the preempt partition.

Mixed systems include a limited number of various GPUs including V100, A100, H100 and GH200.

Up-to-date system information, real-time status, and resource availabilty can be viewed with [OnDemand System Status](https://ondemand-p01.pax.tufts.edu/pun/sys/dashboard/system-status) app.

```{gallery-grid}
---
grid-columns: 1
---
- header: "{fas}`book` Terminology Overview"
  content: "Understand the terminologies used across Tufts HPC cluster documentations"
  link: "../terminology/index.html"

- header: "{fas}`book` HPC Partitions and Limits"
  content: "Partitions on Tufts HPC cluster and their purposes"
  link: "partition.html"

- header: "{fas}`book` GPU Resources"
  content: "Available GPU resources on Tufts HPC cluster and how to request them"
  link: "gpu.html"

```
