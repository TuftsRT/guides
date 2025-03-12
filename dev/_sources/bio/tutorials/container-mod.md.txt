---
tags: bioinformatics
---

# Package installation and module generation with container-mod

Yucheng Zhang, Bioinformatics Engineer, TTS Research Technology
Yucheng.Zhang@tufts.edu

## Introduction

High-Performance Computing (HPC) environments are designed to handle computationally intensive tasks by leveraging powerful hardware, distributed systems, and parallel processing. However, these same features can make package installation more challenging compared to standard desktop or cloud-based computing environments. Users typically lack administrative rights (root access) on HPC systems. This means they cannot use system-wide package managers like apt, yum, or brew to install dependencies.

Tools like Singularity and Apptainer allow users to package their software and dependencies in a container, ensuring portability across HPC systems. However, many users are intimidated by the container syntax.

To simplify the usage of using containers on Tufts HPC, we developed [container-mod](https://github.com/TuftsRT/container-mod). `container-mod` can serve four purposes:

- pull the singularity container image from a public container registry, such as [Docker Hub](https://hub.docker.com/).
- generate the modulefile in lmod format, allowing to use `module load xx` command to use the application.
- generate the wrapper executables for the programs provided by the application.
- [Optional] generate the jupyter kernel for containers (e.g., pytorch and tensorflow) that are able to run on Jupyter Lab/Notebook.

## Syntax

```
container-mod pull|module|exec|pipe [options] URIs
```

## Subcommands

- pull <URI>: Pulls a container image from the provided URI.
- module <URI>: Generates a module file for the container.
- exec <URI>: Creates a wrapper bash script for the container’s programs.
- pipe <URI>: Executes a pipeline that pulls the image, generates a module file, and creates the executable in one step.

## Options

- -d|--dir DIR: Specify the output directory for images, module files, and executables. Defaults to the current directory.
- -f|--force: Force overwrite of existing module files, or executables. Default is to skip existing files.
- -m|--moduledir DIR: Specify the directory that stores module files that can be used as template. Defaults to modulefiles.
- -u|--update: If set, the repository app file will be updated with new version information.
- -p|--personal: Create personal module files in the privatemodules directory (default is no).
- -h|--help: Display this help message and exit.

## Usage

For users, the recommended subcommand is `pipe`, this will run all three subcommands (`pull`, `module`, and `exec`). For options, users need to use `-p` or `--personal`. This will generate the modulefile into users' `$HOME/privatemodules`.
Another commonly used option is `-j` or `--jupyter`. Many users prefer to run python codes within Jupyter lab/notebook on [Tufts Open OnDemand platform](http://ondemand.pax.tufts.edu). However, to enable python codes to run Jupyter lab/notebook, users need to write a jupyter kernel first. These steps are too complex for beginner users. Adding `-j` or `--jupyter`, the jupyter kernel will be created for users. When users start Jupyter lab/notebook, the newly created kernel is ready to use.

## Load the module

```
module load container-mod
container-mod -h
```

## Examples

### vcftools

#### Create the app

```
container-mod pipe -p docker://quay.io/biocontainers/vcftools:0.1.16--h9a82719_5
```

```
Generating executable for vcftools
+-------------------------------------------------+
| To use this module, load the following modules: |
|                                                 |
|     module load use.own                         |
|     module load vcftools/0.1.16                 |
|                                                 |
+-------------------------------------------------+
```

### Use vcftools

```
module load use.own
module load vcftools/0.1.16
```

### pytorch with jupyter support

```
container-mod pipe -p -j docker://tuftsttsrt/pytorch:2.5.1-cuda12.1-cudnn9-runtime-jupyter
```

```
+---------------------------------------------------------------+
| To use this module, load the following modules:               |
|                                                               |
|     module load use.own                                       |
|     module load pytorch/2.5.1-cuda12.1-cudnn9-runtime-jupyter |
|                                                               |
+---------------------------------------------------------------+
Generating Jupyter kernel for pytorch version 2.5.1-cuda12.1-cudnn9-runtime-jupyter
Jupyter kernel created: pytorch-2.5.1-cuda12.1-cudnn9-runtime-jupyter
You can now launch Jupyter Notebook and select the kernel 'pytorch 2.5.1-cuda12.1-cudnn9-runtime-jupyter'
If you'd like to edit the kernel, you can find it at: /cluster/home/tutln02/.local/share/jupyter/kernels/pytorch-2.5.1-cuda12.1-cudnn9-runtime-jupyter
```

### Use pytorch in a jobscript

```
module load use.own
module load pytorch/2.5.1-cuda12.1-cudnn9-runtime-jupyter
```

## Profile mode for the lab/group

**container-mod** allows users to create containerized modules for shared use within their group. The PI or lab/group manager can create a group profile, which can then be used by specifying the `--profile` option on the command line.

### Generate a lab/group profile

Users need to store the profile in `$HOME/container-apps/profiles`.

The following profile `$HOME/container-apps/profiles/rt` for my `rt` group will illustrate the required environment variables.

```
MOD_EXISTING_DIR_DEF="/cluster/tufts/biocontainers/modules"
PUBLIC_IMAGEDIR="/cluster/tufts/rt/shared/container-modules/images"
PUBLIC_EXECUTABLE_DIR="/cluster/tufts/rt/shared/container-modules/tools"
IMG_OUTDIR=$PUBLIC_IMAGEDIR
MOD_OUTDIR="/cluster/tufts/rt/shared/container-modules/modules"
EXEC_OUTDIR=$PUBLIC_EXECUTABLE_DIR
```

- **MOD_EXISTING_DIR_DEF**: container-mod can reuse existing modulefiles. This variable defines the directory where those modulefiles are located. This can be left blank.
- **PUBLIC_IMAGEDIR**: The path to the directory containing singularity images.
- **PUBLIC_EXECUTABLE_DIR**: Directory containing bash wrappers for the target program's commands.
- **IMG_OUTDIR**: This path specifies the directory where singularity images are downloaded. In most cases, this is the same as the directory specified by the **\$PUBLIC_IMAGEDIR** environment variable.
- **MOD_OUTDIR**: The directory in which modulefiles will be stored.
- **EXEC_OUTDIR**: Bash wrappers for the target program's commands are stored here, usually same with **\$PUBLIC_EXECUTABLE_DIR**.

### List available profiles

```
$ container-mod -l   # or container-mod --list
rt  (Personal Profile)
biocontainers
ngc
```

### Run container-mod with lab/group profile

```
$ container-mod pipe --profile rt docker://staphb/samtools:1.21
```

After the workflow is completed, you can see the following message:

```
+----------------------------------------------------------------------------------------------------+
| To use this module, load the following modules:                                                    |
|                                                                                                    |
|     module use /cluster/tufts/rt/shared/container-modules/modules                                  |
|     module load samtools/1.21                                                                      |
|                                                                                                    |
| The modulefile is located at: /cluster/tufts/rt/shared/container-modules/modules/samtools/1.21.lua |
+----------------------------------------------------------------------------------------------------+
```

### Run the lab/group shared module

```
$ module use /cluster/tufts/rt/shared/container-modules/modules
$ module load samtools/1.21
$ samtools view -S -b input.sam > output.bam
$ samtools sort output.bam -o output.sorted.bam
$ samtools index output.sorted.bam
$ samtools flagstat output.sorted.bam
```
