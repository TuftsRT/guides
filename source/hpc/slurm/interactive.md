# Interactive Sessions

An Interactive job is one that a user starts and then continues to directly interact with it during the duration of its execution. This could be any shell application, Jupyter notebooks, X11 graphical application, etc. Interactive jobs are typically submitted via the Slurm salloc or srun commands or the OnDemand web interface.

`$ srun [options] --pty [command]`

- Command

  - `command` to run an application, given the module is already loaded.
  - `bash` to request a bash shell

- Options
  --pty

  - Partition `-p`
    - Default batch
  - Time `-t` or `--time=`
    - Default 15 minutes
    - Format: DD-HH:MM:SS
  - Number of CPU cores `-n`
    - Default 1
    - Common 2 or 4
  - Memory `--mem=`
    - Default 2GB
    - Common 4GB per CPU requested
  - GPU `--gres=`
    - Default none
    - Common --gres=gpu:1
  - X Window `--x11=first`
    - Default none
  - QOS `--qos=`
    - Default normal
    - Common interactive

## CPU Job

Starting an interactive session of bash shell on preempt partition with 2 CPU cores and 2GB of RAM, with X11 forwarding for 1 day, 2 hours, and 30 minutes (use `exit` to end session and release resources).

```bash
[tutln01@login-p01 ~]$ srun -p preempt -t 1-2:30:00 -n 2 --mem=2g --x11=first --pty bash
srun: job 296794 queued and waiting for resources
srun: job 296794 has been allocated resources
[tutln01@pax00# ~]$
```

Note: If you are requesting X11 forwarding with `srun`, `-XC` or`-YC` or `-XYC` must be used upon login with `ssh`.

## GPU Job

Starting an interactive session of bash shell on preempt partition with 2 CPU cores and 4GB of RAM, with 1 A100 GPU with 40GB of device memory for 1 day, 2 hours, and 30 minutes (use `exit` to end session and release resources).

```bash
[tutln01@login-p01 ~]$ srun -p preempt -t 1-2:30:00 -n 2 --mem=4g --gres=gpu:a100:1 --constraint="a100-40G" --pty bash
```

Once your resource is allocated on a compute node, use `nvidia-smi` to check GPU info.
