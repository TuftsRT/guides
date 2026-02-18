# Advanced Jupyter

```{warning}
This tutorial covers an advanced method to install onces own copy of Jupyer.  Most user should utilized our existing Jupyter install via OnDemand.
```

This tutorial is particualrily useful if you need to run on a system that we do not have compatible versions of conda or jupyter in modules for. For example the arm GH200 test system.

## Jupyter

### Install

- These instructions require conda. You can use our conda install available via modules, or the custom conda instructions below.

#. conda create -n jupyter python=3
#. conda activate jupyter
#. conda install -c conda-forge jupyterlab

### Workflow

This is an example work flow to access a jupyter notebook from the browser of your computer.

From a login node request a compute node be allocated

[utln01@login-p01 ~]\$ `salloc -N 1 -p gpu --gres=gpu:1 -t 15`

Once an node is allocated by Slurm you will be connected to this node automatically. You will see your shell prompt change to **[utln01@pax00# ~]\$** .

On this compute node you cam run the commands below to start jupyter

```
conda activate jupyter
jupyter notebook --no-browser --ip=0.0.0.0
```

Connect a new SSH session from your client with a port forward to the allocated node
ssh -L 8888:nodename##:8888 utln@login-prod.pax.tufts.edu

On your client computer use a web browser to navigate to http://127.0.0.1:8888/?token=#########################################

## Conda

This explains how to install your own local version of conda.

### Installing conda on x86

- wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
- bash Miniconda3-latest-Linux-x86_64.sh -p ~/barn/miniconda3x86
- You will need to read and accept the license terms.
- If this is the first conda install inside your project you can let conda auto initialize your environment.
  -- When prompted to run conda init, answer yes.
  -- source ~/.bashrc
- If this is a 2nd install of conda for a difference architecture do not run conda init. The installer will display a command to run to manually setup conda each time you login.

### Installing conda on Arm

- wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh
- bash Miniconda3-latest-Linux-aarch64.sh
- You will need to read and accept the license terms.
- If this is the first conda install inside your project you can let conda auto initialize your environment.
  -- When prompted to run conda init, answer yes.

### Common situations

#### Manually activating conda, or switching conda between x86/arm installs

This command will activate the conda install at the path provided in your current shell `eval "$(~/miniconda3x86/bin/conda shell.bash hook)"`

#### Activating a conda environment in a script or Slurm submission file

You can a activate conda environment using the following commands. Replace the path in source with the path for you conda install

```
source ~/miniconda3x86/etc/profile.d/conda.sh
conda activate environment-name
```

#### Running job using sbatch

Sometimes it makes sense to run Jupyter as a "non interactive" job using sbatch instead of salloc. This can be useful for long running notebooks or other situations where you want the job to keep running if you are disconnected. Even though the job is "non interactive" you can still connect to it the same way and use it interactively. There are 3 main differences when using this method.

**Main Differences**

1. Use sbatch and a job submission script
1. Locate your allocated node
1. Get Jupyter port and URL from output log

##### Use sbatch and a job submission script

To start your job create a file named jupyter_session.sh

```
#!/bin/bash
#SBATCH -p PartitionName  # gpu, batch or preempt
#SBATCH -t 04:00:00 #4 Hours
#SBATCH -c 4 #4 CPUs
#SBATCH --gres=gpu:1 # 1 GPU, if needed
#SBATCH --job-name=jupyter
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --error=%x-%J-%u.err
#SBATCH --output=%x-%J-%u.out

module purge
module load conda

conda activate jupyter
jupyter notebook --no-browser --ip=0.0.0.0
```

##### Locate your allocated node

The squeue command will show you your jobs, and the nodes they are allocated to run on. Use `squeue --me` and locate the hostname (pax###) of the node that shows your jobs as running.

`squeue --me`

##### Establish SSH port forward

Setup the SSH port forward
