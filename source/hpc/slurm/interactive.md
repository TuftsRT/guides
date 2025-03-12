# Interactive Sessions

**Particularly good for debugging and working with software GUI.**

We strongly recommend user check available resources using [hpctools - Tufts HPC Helper Tool](../examples/hpctools.md) before starting interactive jobs.

`$ srun [options] --pty [command]`

- Command

  - command to run an application, given the module is already loaded.
  - `bash` for a bash shell

- Options

  - Pseudo terminal `--pty`
  - Partition `-p`
    - Default batch if not specified
  - Time `-t` or `--time=`
    - **Default 15 minutes** if not specified on non-interactive partition
    - Format: DD-HH:MM:SS
  - Number of CPU cores `-n`
    - Default 1 if not specified
  - Memory `--mem=`
    - Default 2GB if not specified
  - GPU `--gres=`
    - Default none
  - X Window `--x11=first`
    - Default none
  - Specify node name `-w` or `--nodelist=`
    - Job will only be placed on the nodes listed.
    - Optional.

## CPU Job

**Starting an interactive session of bash shell on preempt partition with 2 CPU cores and 2GB of RAM, with X11 forwarding for 1 day, 2 hours, and 30 minutes (use `exit` to end session and release resources).**

```bash
[tutln01@login-prod-01 ~]$ srun -p preempt -t 1-2:30:00 -n 2 --mem=2g --x11=first --pty bash
srun: job 296794 queued and waiting for resources
srun: job 296794 has been allocated resources
[tutln01@cc1gpu001 ~]$
```

Note: If you are requesting X11 forwarding with `srun`, `-XC` or`-YC` or `-XYC` must be used upon login with `ssh`.

## GPU Job

**Starting an interactive session of bash shell on preempt partition with 2 CPU cores and 4GB of RAM, with 1 A100 GPU with 40GB of device memory for 1 day, 2 hours, and 30 minutes (use `exit` to end session and release resources).**

```bash
[tutln01@login-prod-01 ~]$ srun -p preempt -t 1-2:30:00 -n 2 --mem=4g --gres=gpu:a100:1 --constraint="a100-40G" --pty bash
```

Once your resource is allocated on a compute node, use `nvidia-smi` to check GPU info.
