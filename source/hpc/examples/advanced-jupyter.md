# Advanced Jupyter 

```{warning}
This tutorial covers an advanced method to install onces own copy of Jupyer.  Most user should utilized our existing Jupyter install via OnDemand.
```

This tutorial is particualrily usefull if you need to run on a system that we do not have compatible versions of conda or jupyter in modules for.  For example the arm GH200 test system.

## Jupyter

### Install
* These instructions require conda.  You can use conda via modules, or the custom conda instructions below.

#. conda create -n jupyter python=3
#. conda activate jupyter
#. conda install -c conda-forge jupyterlab

### Workflow
This is an example work flow to access a jupyter notebook from the browser of your computer.

Request a node be allocated

`salloc -N 1 --gres=gpu:1 -t 15`

See which node was allocated

`squeue`

SSH to the node and start jupyter notebook
```
ssh nodename##
conda activate jupyter
jupyter notebook --no-browser --ip=0.0.0.0
```

Connect a new SSH session from your client with a port forward to the allocated node
ssh -L 8888:nodename##:8888 username@login.pax.tufts.edu

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
