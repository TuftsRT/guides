# Advanced Jupyter

```{warning}
This tutorial covers an advanced method to install one's own copy of Jupyer.  Most users should utilized our existing Jupyter install via OnDemand.
```

## Who is this tutorial for

This tutorial is particularly useful if you need to run on a system that we do not have compatible versions of conda or Jupyter, such as ARM-based systems like the GH200, or if you want to access a running Jupyter session from your local machine.

### Port forwarding

This also covers port forwarding, a process that allows users to access Jupyter running on the Tufts HPC Cluster from their local computer. This is especially useful for connecting to the LLM Notebooks from your local computer, to allow for running Ollama models on the HPC from your local machine's familiar environment. For more details on the AI Tools on the cluster, see here: [https://go.tufts.edu/AITools](../../ai/index).

## Jupyter

## Install

These instructions require conda. You can use our [conda install available via modules](../application/10-condaenv), or the [custom conda instructions below](#conda).

```
conda create -n jupyter python=3
conda activate jupyter
conda install -c conda-forge jupyterlab
```

## Workflow

This is an example work flow to access a Jupyter notebook from the web browser of your local computer.

Using SSH on your local computer, connect to a login node and request a compute node be allocated

[utln01@login-p01 ~]\$ ```salloc -N 1 -p gpu --gres=gpu:1 -t 15```

Once an node is allocated by Slurm you will be connected to this node automatically. You will see your shell prompt change to **[utln01@pax00# ~]\$** .

On this compute node you cam run the commands below to start Jupyter

```
conda activate jupyter
jupyter notebook --no-browser --ip=0.0.0.0
```

As Jupyter starts up you will see several message printed to your SSH console.  From this you need 3 important pieces of information.
1. The network port that Jupyter started on.  This is a 4 digit number, typically `8888`
2. The local host (127.0.0.1) URL with authentication token Jupyter generated.  It will look similar to `http://127.0.0.1:8888/?token=#########################################`.
3. The nodename of the compute node you are running Jupyter on. This will be pax with 3 numbers after it, similar to pax000.  It will be visible in the message, and is also the same as the name shown in the shell prompt.

In a new window, connect a 2nd SSH session from your client with a port forward to the allocated node
```
ssh -L 8888:pax###:8888 utln@login-prod.pax.tufts.edu
```

On your client computer use your favorite web browser to navigate to http://127.0.0.1:8888/?token=######################################### replacing this example url with the similar localhost (127.0.0.1) version output by Jupyter when it started.

## Conda

This explains how to install your own local version of conda.

### Installing conda on x86

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -p ~/miniconda3x86
```
- You will need to read and accept the license terms.
- If this is the first conda install inside your account you can let conda auto initialize your environment.
  -- When prompted to run conda init, answer yes.
  -- source ~/.bashrc
- If this is a 2nd install of conda for a difference architecture do not run conda init. The installer will display a command to run to manually setup conda each time you login.  Keep this for future reference.

### Installing conda on Arm

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh
bash Miniconda3-latest-Linux-aarch64.sh  -p ~/miniconda3arm
```
- You will need to read and accept the license terms.
- If this is the first conda install inside your account you can let conda auto initialize your environment.
  -- When prompted to run conda init, answer yes.

## Common situations

### Manually activating conda, or switching conda between x86/arm installs

This command will activate the conda install at the path provided in your current shell.
 ```eval "$(~/miniconda3x86/bin/conda shell.bash hook)"```

### Activating a conda environment in a script or Slurm submission file

You can activate a conda environment using the following commands. Replace the path passed to the source command with the path for you conda install

```
source ~/miniconda3x86/etc/profile.d/conda.sh
conda activate environment-name
```

### Different port number than 8888

It is not uncommon for multiple instances of Jupyter to be running on the same compute node. If this is the case, the **default** port number of 8888 will be a different number, you will need to use this alternate number instead.

You can specify a port number using the `--port` option when calling `jupyter lab` or `jupyter notebook`, for example `--port 6767`. However, Jupyter will fail to start if that number you manually specify is already in use.

To search for available ports, you can use the command:
```ss -tuln```

We recommend you use ports around 8000 for this. You can check if a given port is open, such as 8080, by passing the result to `grep`:

```ss -tuln | grep -w ':8080'```

This will return a service it is exists on the port, and nothing if the port is available.

### Starting Jupyter using sbatch

Sometimes it makes sense to run Jupyter as a "non interactive" job using sbatch instead of salloc. This can be useful for long running notebooks or other situations where you want the job to keep running if you are disconnected. Even though the job is "non interactive" you can still connect to it the same way and use it interactively. There are 3 main differences when using this method.

**Main Differences**

1. Use sbatch and a job submission script
1. Locate your allocated node
1. Get Jupyter port and URL from output log

#### Use sbatch and a job submission script

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
module load anaconda/2025.06.0

conda activate jupyter
jupyter lab --no-browser --ip=0.0.0.0
```

#### Locate your allocated node

The squeue command will show you your jobs, and the nodes they are allocated to run on. Use `squeue --me` and locate the hostname (pax###) of the node that shows your jobs as running.

```squeue --me```

```
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
            247071       gpu  jupyter utln01  R       0:06      1 pax0##
```

#### Locating connection information

When using sbatch the Jupyter output will be located in the slurm error log, in this case a file in the same folder as your submission script in a file named similar to jupyter-jobid#-utln01.err . This file will contain familiar output that includes a url similar to the one shown below.

`http://127.0.0.1:8888/lab?token=234234asadf959invalid96223453dsa4de`

Setup the SSH port forward using the node name we found previously and the port number shown in this url found in the error log.

```
ssh -L 8888:pax###:8888 utln@login-prod.pax.tufts.edu
```

On your client computer use a web browser to navigate to http://127.0.0.1:8888/lab?token=234234asadf959invalid96223453dsa4de
