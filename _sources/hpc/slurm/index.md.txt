# Slurm Job Scheduler

The HPC cluster uses the Slurm Job Scheduler to assign users jobs to compute nodes. Jobs are allocated based on the
requested resources and the submitting users priority. We use the FairShare algorithm which adjusts priority to
balance usage across of our users.

Command Quick Reference

- **squeue** lists your jobs in the queue
- **sinfo** lists the state of all computers in the HPC cluster
- **sbatch** submits batch jobs
- **sprio** Displays the priorities of pending jobs in the queue
- **scancel** can be used to cancel jobs
- **salloc** allocates a compute node for interactive use
- **sacct** display historical report data for jobs
- **seff** displays job CPU and Memory efficiency

```{gallery-grid}
---
grid-columns: 1
---
- header: "{fas}`book` Partitions"
  content: "How to select where your jobs and the resources available to them."
  link: "../compute/partition.html"

- header: "{fas}`book` Interactive Sessions"
  content: "How to run interactive jobs with cluster resources"
  link: "interactive.html"

- header: "{fas}`book` Batch jobs"
  content: "How to submit batch jobs to the cluster and sample slurm batch job scripts"
  link: "batchjob.html"

- header: "{fas}`book` Job Monitoring and Management Commands"
  content: "How to monitor and manage your active jobs"
  link: "monitor.html"

- header: "{fas}`book` Job Resource Utilization"
  content: "How to check resource utilization of completed jobs: how much memory did the completed job used? how many cores did the completed job utilized?"
  link: "utilization.html"


```
