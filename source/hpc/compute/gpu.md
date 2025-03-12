# GPU Resources

GPUs (Graphics Processing Units) are essential for tasks requiring high parallelism, such as deep learning and simulations. In the context of HPC these devices no longer perform any graphics function, but the "graphics" name has stuck vs the more apt accelerator term. To lean more about accelerated computing visit https://blogs.nvidia.com/blog/what-is-accelerated-computing/.

NVIDIA GPUs are available in `gpu` and `preempt` public partitions as well as some labs private partitions. When scheduling batch or interactive jobs that need GPU you add the `--gres` option to your command. This allows you to select the type and quantity of GPUs needed.

- The simplest request is `--gres=gpu:1`, e.g.

  `$ srun -p preempt -n 2 --mem=4g --gres=gpu:1 -t 1:00:00 --pty bash`

  - `--gres` : Generic Resource
  - `gpu:1` : requesting one GPU (any GPU architecture available in the request partition)

- To request a specific GPU architecture you can add it to the gres `--gres=gpu:t4:1`, e.g.

  `$ srun -p preempt -n 2 --mem=4g --gres=gpu:t4:1 -t 1:00:00 --pty bash`

  - `--gres` : Generic Resource
  - `gpu:t4:1` : requesting one T4 GPU.

- You can request more than one type (but not all) of GPUs with `constraint`, e.g.

  `$ srun -p preempt -n 2 --mem=4g --gres=gpu:1 --constraint="t4|p100|v100" -t 1:00:00 --pty bash`

  - `--constraint` : set constraints on the resources allocated for the task.
  - `t4|p100|v100` : indicates that the task can use a GPU of type t4, p100, or v100, where the | symbol means "or".

```{warning}
- **DO NOT** manually set `CUDA_VISIBLE_DEVICES`.

- Users can only see GPU devices that are assigned to them by using the `$ nvidia-smi` command.

- When submitting batch jobs, it's recommended adding this command `$ nvidia-smi` in your slurm job submission script to
include the output in the log for troubleshooting purposes.
```

## Available GPUs

| GPU Model   | Memory     | Max per node | Partition Available in | Constraints                     | Notes                    |
| ----------- | ---------- | ------------ | ---------------------- | ------------------------------- | ------------------------ |
| a100        | 40G or 80G | 8            | gpu, preempt           | a100 <br/>a100-40G<br/>a100-80G |                          |
| p100        | 16GB       | 6            | gpu, preempt           | p100                            |                          |
| v100        | 16GB       | 4            | preempt                | v100                            |                          |
| t4          | 16GB       | 4            | preempt                | t4                              | Max CUDA version is 10.2 |
| rtx_6000    | 24GB       | 8            | preempt                | rtx_6000                        |                          |
| rtx_a6000   | 48GB       | 8            | preempt                | rtx_a6000                       |                          |
| rtx_6000ada | 48GB       | 4            | preempt                | rtx_6000ada                     |                          |
| l40s        | 48GB       | 4            | preempt                | l40s                            |                          |
| rtx_a5000   | 24GB       | 8            | preempt                | rtx_a5000                       |                          |
| l40         | 48GB       | 4            | preempt                | l40                             |                          |
| h100        | 80GB       | 3            | preempt                | h100                            |                          |

All GPU cards drivers except t4 support CUDA 12.2
