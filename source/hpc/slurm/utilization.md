# Job Resource Utilization

Making sure you are making efficient use of the HPC resources you request is very important.  Any HPC job will only use a portion of the resources you request, so it is important to make sure you are not requesting alot more than you need. Reviewing your completed jobs helps you make better decisions on how much resources to request for future jobs.

Way to view job effiency after completion
- Slurm job completion emails
- Open OnDemand Active Jobs Dashboard
- Slurm seff command

Ways to view job effiency during execution
- top command
- nvidia-smi command

## Completed Jobs

### Slurm job completion emails
If you use Open OnDemand to run jobs, or request Slurm emails on job completion these message include **Job Efficiency Metrics** in the message.  An email attachment includes additional details.

### Active Jobs Dashboard
The Open OnDemand server provides an ["Active Jobs" dashboard](https://ondemand-prod.pax.tufts.edu/pun/sys/dashboard/activejobs) that shows your job history.  It can be accessed at under the **Jobs** -> **Active Jobs** menu item.  Clicking on each completed job in the **Job History** section will show additional details about it.  Completed jobs will include efficiency data under the **Resources and I/O** tab.

### seff 

You can check your **completed** jobs using the slurm seff command.  The would be job that no longer show as running or pending in the queue.  This command can be run on any hpc node. Knowing your slurm JOBID is required.

Display job CPU and memory usage:

`$ seff JOBID`

```bash
[tutln01@login-p01 ~]$ seff 296794
Job ID: 296794
Cluster: pax
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

## Running Jobs

Sometimes it is neccessary to review the resource utilization while a job is running.  This is typically done by connected to the node your jobs is currently running on via SSH, and running additional commands.

### top
The top command will list all the processes running on the compute node and the memory and CPU being used.  To see the processes for your user run `top -u $USER`.

### nvidia-smi
When using GPU nodes the `nvidia-smi` command is useful for spot checking GPU compute and memory usage. You will only see the GPUs allocated to your user.

```
[utln@pax00# ~]$ nvidia-smi
Wed Feb 18 22:52:34 2026       
+-----------------------------------------------------------------------------------------+
| NVIDIA-SMI 575.57.08              Driver Version: 575.57.08      CUDA Version: 12.9     |
|-----------------------------------------+------------------------+----------------------+
| GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
| Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
|                                         |                        |               MIG M. |
|=========================================+========================+======================|
|   0  NVIDIA H200                    On  |   00000000:19:00.0 Off |                    0 |
| N/A   44C    P0            263W /  700W |   18193MiB / 143771MiB |     54%      Default |
|                                         |                        |             Disabled |
+-----------------------------------------+------------------------+----------------------+
|   1  NVIDIA H200                    On  |   00000000:3B:00.0 Off |                    0 |
| N/A   35C    P0             78W /  700W |       0MiB / 143771MiB |      0%      Default |
|                                         |                        |             Disabled |
+-----------------------------------------+------------------------+----------------------+
|   2  NVIDIA H200                    On  |   00000000:4C:00.0 Off |                    0 |
| N/A   43C    P0            252W /  700W |   18711MiB / 143771MiB |     53%      Default |
|                                         |                        |             Disabled |
+-----------------------------------------+------------------------+----------------------+
|   3  NVIDIA H200                    On  |   00000000:5D:00.0 Off |                    0 |
| N/A   45C    P0            257W /  700W |   18711MiB / 143771MiB |     53%      Default |
|                                         |                        |             Disabled |
+-----------------------------------------+------------------------+----------------------+
|   4  NVIDIA H200                    On  |   00000000:9B:00.0 Off |                    0 |
| N/A   40C    P0            237W /  700W |   18711MiB / 143771MiB |     47%      Default |
|                                         |                        |             Disabled |
+-----------------------------------------+------------------------+----------------------+
|   5  NVIDIA H200                    On  |   00000000:BB:00.0 Off |                    0 |
| N/A   34C    P0             79W /  700W |       0MiB / 143771MiB |      0%      Default |
|                                         |                        |             Disabled |
+-----------------------------------------+------------------------+----------------------+
|   6  NVIDIA H200                    On  |   00000000:CB:00.0 Off |                    0 |
| N/A   42C    P0            241W /  700W |   18711MiB / 143771MiB |     48%      Default |
|                                         |                        |             Disabled |
+-----------------------------------------+------------------------+----------------------+
|   7  NVIDIA H200                    On  |   00000000:DB:00.0 Off |                    0 |
| N/A   34C    P0             80W /  700W |       0MiB / 143771MiB |      0%      Default |
|                                         |                        |             Disabled |
+-----------------------------------------+------------------------+----------------------+
                                                                                         
+-----------------------------------------------------------------------------------------+
| Processes:                                                                              |
|  GPU   GI   CI              PID   Type   Process name                        GPU Memory |
|        ID   ID                                                               Usage      |
|=========================================================================================|
|    0   N/A  N/A         1213801      C   python                                18184MiB |
|    2   N/A  N/A         1213800      C   python                                18702MiB |
|    3   N/A  N/A         1213803      C   python                                18702MiB |
|    4   N/A  N/A         1059301      C   python                                18702MiB |
|    6   N/A  N/A         1213802      C   python                                18702MiB |
+-----------------------------------------------------------------------------------------+
```