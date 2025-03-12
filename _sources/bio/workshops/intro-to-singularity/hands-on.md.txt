# Running Containerized Applications on Tufts HPC: A Guide to Pulling and Executing

## Start an interactive job

Users cannot run jobs on login nodes in an HPC system because login nodes are designed for lightweight tasks such as:

- Editing scripts
- Submitting jobs to the slurm scheduler
- Managing files/folders

To pull a container or run a container, start an interactive session on a compute node:

```
$ srun -N1 -n2 -t1:00:00 -p batch --pty bash ## You will start an interactive session on batch partition with 2 cores and 2-hour walltime.
```

## Load modules

It is highly recommended to run `module purge` before loading `singularity`.

```
$ module purge
$ module load singularity/3.8.4
$ module list

Currently Loaded Modules:
  1) squashfs/4.4   2) singularity/3.8.4
```

Sure, here's an improved version:

`squashfs/4.4` is a necessary dependency for `singularity`. When you load `singularity/3.8.4`, `squashfs/4.4` will be automatically loaded as well.

## singularity pull

### Syntax

Download or pull a container from a given URI.

```
$ singularity pull [output file] <URI>
```

### Samtools

BioContainers is integrated with bioconda. You can find almost all bioinformatics applications from BioContainer's [Package Index](https://bioconda.github.io/conda-package_index.html)

[samtools](https://github.com/samtools/samtools) is an bioinformaitcs software for dealing with SAM, BAM and CRAM files. Here I will show how to pull its latest version (1.21.0) from BioContainers.

Here is the samtools's [bioconda page](https://bioconda.github.io/recipes/samtools/README.html#package-package%20'samtools').

There are different tags for each version and multiple containers/tags even exist for the same application version. It is common practice to select the ones that were last modified.

To pull the image from BioContainers, we can use the following command:

```
## Default
$ singularity pull docker://quay.io/biocontainers/samtools:1.21--h96c455f_1

## To give a customerized output name
$ singularity pull samtools_1.21.sif docker://quay.io/biocontainers/samtools:1.21--h96c455f_1
```

### PyTorch

[PyTorch](https://pytorch.org) is a powerful open-source machine learning framework based on the Python programming language and the Torch library. It's widely used for deep learning, a type of machine learning that builds complex models like artificial neural networks for tasks like image recognition, natural language processing, and more.

If you want to install PyTorch, it can be complex. You have to ensure different libraries (e.g., CUDA) and dependencies are compatible. If you use container, you can easily an container rom Docker Hub. Here we want to show you how to pull the pytorch container under Tufts Research Technology account [tuftsttsrt](https://hub.docker.com/r/tuftsttsrt/pytorch/tags):

```
$ singularity pull docker://tuftsttsrt/pytorch:2.5.1-cuda12.1-cudnn9-runtime-jupyter
```

## singularity exec

Singularity has a `run` subcommand. This can be used to run the user-defined default command within a container. However, it is recommended to use another subcommand `exec` instead of `run`. If you're interested in `run`, please check the [singularity user guide](https://docs.sylabs.io/guides/3.8/user-guide/cli/singularity_run.html#singularity-run).

### Syntax

Download or build a container from a given URI.

```
singularity exec [options] image command
```

### samtools

```
$ singularity exec samtools_1.21--h96c455f_1.sif samtools
```

### pytorch

```
$ singularity exec pytorch_2.5.1-cuda12.1-cudnn9-runtime-jupyter.sif python
Python 3.11.10 | packaged by conda-forge | (main, Oct 16 2024, 01:27:36) [GCC 13.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import torch
>>> print(torch.__version__)
2.5.1+cu121
```

## container-mod

### samtools

```
$ module load container-mod
$ container-mod pipe docker://quay.io/biocontainers/samtools:1.21--h96c455f_1
$ module load use.own
$ module load samtools/1.21
```

### pytorch

```
$ module load container-mod
$ container-mod pipe -j docker://tuftsttsrt/pytorch:2.5.1-cuda12.1-cudnn9-runtime-jupyter
$ module load use.own
$ module load pytorch/2.5.1-cuda12.1-cudnn9-runtime-jupyter
```

- `-j` is used when the container supports jupyter. With `-j`, the jupyter kernel will be created, allowing you to run `pytorch` on **Jupyter Lab/Notebook** on Tufts [Open OnDemand](http://ondemand.pax.tufts.edu).
