# Conda Environments

Conda is a software package manager similar to PIP that allows users to easily install software not already available on the cluster. It is designed to install software into user space and does not require administrative permissions. This makes it ideal for use on shared HPC systems.

In most cases `conda install package-name` is a direct replacement for PIP in instructions from the internet you might be following. You can search [Anaconda Cloud](https://anaconda.org/search) for packages available for install. Instructions on specifying an exact version can be found at [Managing Packages](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-pkgs.html).

## Allocate Resources

1. [Start an interactive session](../slurm/interactive) on the cluster on a compute node (execute every time you need to execute/run any programs on the cluster)

   1. Determine if you need GUI forwarding
   1. Determine if you need GPU access

1. **Load modules**

   1. Check available conda modules

      `$ module av miniforge`

      `$ module av anaconda`

      `$ module av miniconda`

   1. Load module

      `$ module load miniforge/24.11.2-py312` **Recommended**

      `$ module load anaconda/2021.05`

      or any newer versions of anaconda shown in the `module av` output

      NOTE: DO NOT use `anaconda/3` and `anaconda/2`, they are very old versions.

   1. Load other modules needed (such as `$ module load cuda/12.2`)

## Configure your conda

**_NOTE: in most cases, the steps in this session only needs to be executed ONCE_**

Since you have limited amount of storage in your home directory, we do no suggest you install the packages there. As you belong to XXXXlab group on the cluster, please use the group research storage for the purpose.

Create two directories in your group research storage space (one for storing the envs, one for storing the pkgs, for example: condaenv, condapkg)

`$ mkdir /cluster/tufts/XXXXlab/$USER/condaenv/`

`$ mkdir /cluster/tufts/XXXXlab/$USER/condapkg/`

If you haven't used conda before on the cluster, create a file named ".condarc" in your home directory.

Now add the following 4 lines to the `.condarc` file in your home directory (modify according to your real path to the directories):

```
envs_dirs:
  - /cluster/tufts/XXXXlab/$USER/condaenv/
pkgs_dirs:
  - /cluster/tufts/XXXXlab/$USER/condapkg/
```

**OR** you can do so from command line with the following commands

`$ conda config --append envs_dirs /cluster/tufts/XXXXlab/$USER/condaenv/`

`$ conda config --append pkgs_dirs /cluster/tufts/XXXXlab/$USER/condapkg/`

**_(only add 4 lines to the .condarc file OR use these commands, NOT BOTH)_**

Optional: Add channels to your conda config as well (here are only 2, but you can add more you need):

`$ conda config --add channels bioconda`

`$ conda config --add channels conda-forge`

After this, your `.condarc`file should look like this

`$ cat ~/.condarc`

```
envs_dirs:
  - /cluster/tufts/XXXXlab/$USER/condaenv/
pkgs_dirs:
  - /cluster/tufts/XXXXlab/$USER/condapkg/
channels:
  - bioconda
  - conda-forge
  - defaults
```

## Create your conda environment

Now you can create your own conda env

`$ cd /cluster/tufts/XXXXlab/$USER/condaenv/`

`$ conda create -p yourenvname`

or if you have a specific version of python you need to use, e.g. 3.8

`$ conda create -p yourenvname python=3.8` (Recommended!)

Note: you will need to have `python` and `pip` installed inside the env to pip install packages inside the env.

Activate the environment (needs to be executed whenever you need to use the conda env you have created)

`$ source activate yourenvname`

If you are using system installed conda, please DO NOT use `conda activate` to activate your environment

Install `yourpackage` in the conda env

`$ conda install yourpackage`

Or if you have python (comes with pip) installed,

`$ pip install yourpackage`

Or follow the instruction on package website.

Check what's installed in your conda environment

`$ conda list`

Test

`$ python`

`<<< import yourpackage`

`<<<`

When you are done, deactivate the environment

`$ source deactivate`

or

`$ conda deactivate`

:)

## Additional Information for Jupyter Users

### Run conda env as a kernel in Jupyter

If you would like to use JupyterNotebook or JupyterLab from OnDemand, you can follow the instructions below and **run your conda env as a kernel in Jupyter**.

- Make sure with python 3.7+ and make sure you load cluster's anaconda module (this only works with py3.7+)

- Activate your conda env from terminal.

- Install ipykernel with

  `$ pip install ipykernel` (this assumes you installed python and pip in your env, otherwise, use "--user" flag)

- Add your env to jupyter with

  `$ python -m ipykernel install --user --name=myenvname` (use the name of your environment in the place of "myenvname")

- Restart Jupyter from OnDemand

- :)
