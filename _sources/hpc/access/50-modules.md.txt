# Modules

## What are modules?

- A tool that **simplify** the use of different software (versions) in a precise and controlled manner by initializing and configuring your shell environment without manual intervention. They allow users easily modify their shell environment during the session with **modulefiles**.

- Modulefiles

  - Each modulefile contains the **information** needed to configure the shell environment for an application. (PATH, LD_LIBRARY_PATH, CPATH, etc).

  - By loading a modulefile, the system automatically adjusts these variables so that the chosen application is ready to use.

  - Modules are useful in managing **different versions** of applications.

  - Modules can also be bundled into meta-modules that will load an entire **set of different applications**.

## Why Use Modules?

- Using modules simplifies managing your computational environment, making it easy to ensure that all necessary software and their correct versions are available for your tasks. This reduces potential conflicts and errors, streamlining your workflow and improving efficiency on the HPC cluster.

## How to use modules?

> **Module Cheat Sheet**
>
> `module av`
>
> `module av <software>`
>
> `module list`
>
> `module load <software>`
>
> `module unload <software>`
>
> `module swap <loaded_software> <new_software>`
>
> `module purge`

- To **check available modules** installed on the cluster:

  ```
   [tutln01@login-prod-01 ~]$ module av
  ```

- Upon login, environment variable **`PATH`** is set for the system to search executables along these paths:

  ```
  [tutln01@login-prod-01 ~]$ echo $PATH

  /usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/cluster/home/tutln01/bin:/cluster/home/tutln01/.local/bin
  ```

### Example: Using `gcc` Compiler

- Check what versions of gcc compiler is available: load the version I would like to use, and use it:

  ```
  [tutln01@login-prod-01 ~]$ module av gcc

  ----------------------------------------------------------- /opt/shared/Modules/modulefiles-rhel6     ------------------------------------------------------------
  gcc/4.7.0 gcc/4.9.2 gcc/5.3.0 gcc/7.3.0

  -------------------------------------------------------------- /cluster/tufts/hpc/tools/module   ---------------------------------------------------------------
  gcc/8.4.0 gcc/9.3.0 gcc/11.2.0
  ```

- Load the desired version of gcc:

  ```
  [tutln01@login-prod-01 ~]$ module load gcc/7.3.0
  ```

- Use `module list` to **check loaded modules** in current environment:

  ```
  [tutln01@login-prod-01 ~]$ module list

  Currently Loaded Modulefiles:
  1) use.own     2) gcc/7.3.0
  ```

- Verify the `gcc` Version:

  ```
  [tutln01@login-prod-01 ~]$ which gcc
  /opt/shared/gcc/7.3.0/bin/gcc

  [tutln01@login-prod-01 ~]$ echo $PATH
  /opt/shared/gcc/7.3.0/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/cluster/home/tutln01/bin:/cluster/ho  me/tutln01/.local/bin

  [tutln01@login-prod-01 ~]$ gcc --version
  gcc (GCC) 7.3.0
  Copyright (C) 2017 Free Software Foundation, Inc.
  This is free software; see the source for copying conditions.  There is NO
  warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  ```

- **swap a module for another** (doesn't have to be the same software):

  ```
  [tutln01@login-prod-01 ~]$ module swap gcc/7.3.0 gcc/9.3.0
  [tutln01@login-prod-01 ~]$ module list
  Currently Loaded Modulefiles:
  1) use.own     2) gcc/9.3.0
  ```

- **unload loaded modules**:

  ```
  [tutln01@login-prod-01 ~]$ module unload gcc
  [tutln01@login-prod-01 ~]$ echo $PATH
  /usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/cluster/home/tutln01/bin:/cluster/home/tutln01/.local/bin
  ```

- **unload ALL** of the loaded modules:

  ```
  [tutln01@login-prod-01 ~]$ module purge
  ```

​

## Private Modules

Users can create [Private Modules](../examples/privatemodule) for the software installed locally.
