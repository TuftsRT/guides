# Job Resource Utilization

## Check your **finished** jobs

_The jobs you can no longer see in `squeue` command output._

**Querying finished jobs helps users make better decisions on how much resources to request for future jobs.**

JOBID is required.

Display job CPU and memory usage:

`$ seff JOBID`

```bash
[tutln01@login-prod-01 ~]$ seff 296794
Job ID: 296794
Cluster: pax
Use of uninitialized value $user in concatenation (.) or string at /usr/bin/seff line 154, <DATA> line 602.
User/Group: /tutln01
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 2
CPU Utilized: 00:00:00
CPU Efficiency: 0.00% of 00:22:12 core-walltime
Job Wall-clock time: 00:11:06
Memory Utilized: 1.16 MB (estimated maximum)
Memory Efficiency: 0.06% of 2.00 GB (2.00 GB/node)
```

Display job detailed accounting data:

`$ sacct --format=partition,state,time,start,end,elapsed,MaxRss,ReqMem,MaxVMSize,nnodes,ncpus,nodelist -j JOBID`

```bash
[tutln01@login-prod-01 ~]$ sacct --format=partition,state,time,start,end,elapsed,MaxRss,ReqMem,MaxVMSize,nnodes,ncpus,nodelist -j  296794
 Partition      State  Timelimit               Start                 End    Elapsed     MaxRSS     ReqMem  MaxVMSize   NNodes      NCPUS        NodeList
---------- ---------- ---------- ------------------- ------------------- ---------- ---------- ---------- ---------- -------- ---------- ---------------
   preempt  COMPLETED 1-02:30:00 2021-03-22T22:18:55 2021-03-22T22:30:01   00:11:06                   2Gn                   1          2       cc1gpu001
           OUT_OF_ME+            2021-03-22T22:18:55 2021-03-22T22:30:01   00:11:06         8K        2Gn    135100K        1          2       cc1gpu001
            COMPLETED            2021-03-22T22:18:56 2021-03-22T22:30:01   00:11:05       592K        2Gn    351672K        1          2       cc1gpu001
```

NOTE: there are more format options, see [sacct](https://slurm.schedmd.com/sacct.html)

```{Hint}
If you don't know your job IDs or can't find them.  Utilize [hpctools - Tufts HPC Helper Tool](../examples/hpctools.md) on the HPC cluster to find all of your jobs since a certain date.
```
