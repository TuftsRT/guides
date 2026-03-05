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
> `module spider <software>`
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
   $ module av
  ```

  or

  ```
  $ module spider
  ```

### `module av` vs. `module spider`

| Feature             | `module av`                                                               | `module spider`                                                                  |
| :------------------ | :------------------------------------------------------------------------ | :------------------------------------------------------------------------------- |
| **Search Scope**    | Only shows modules currently available to load based on your environment. | Searches the entire software hierarchy, including hidden or nested modules.      |
| **Hierarchy Aware** | Yes; it only shows what matches your current compiler/MPI stack.          | No; it finds software regardless of whether its dependencies are loaded.         |
| **Common Use Case** | Checking what you can use _right now_.                                    | Finding _if_ a piece of software exists and how to load it.                      |
| **Output Detail**   | Lists available versions briefly.                                         | Provides detailed instructions on required parent modules for specific versions. |

### Example: Using `gcc` Compiler

- Check what versions of gcc compiler is available: load the version I would like to use, and use it:

  ```
  $ module avail gcc
  ------ /cluster/tufts/apps/modules/9/x86_64/Core -------
    gcc/12.4.0    gcc/15.1.0 (D)

  Where:
  D:  Default Module
  ```

- Load the desired version of gcc:

  ```
  $ module load gcc/12.4.0
  ```

- Use `module list` to **check loaded modules** in current environment:

  ```
  $ module list

  Currently Loaded Modules:
  1) gmp/6.3.0-gt7twni (H)   2) gcc/12.4.0
  ```

- Verify the `gcc` Version:

  ```
  $ which gcc
  /cluster/tufts/apps/spack/9/x86_64/apps/linux-broadwell/gcc-12.4.0-c4qsn6yjbnde4mdnr6wlplee454jnden/bin/gcc

  $ gcc --version
  gcc (Spack GCC) 12.4.0
  Copyright (C) 2022 Free Software Foundation, Inc.
  This is free software; see the source for copying conditions.  There is NO
  warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  ```

- **switch a module for another** (doesn't have to be the same software):

  ```
  $ module switch gcc/12.4.0 gcc/15.1.0
  $ module list
  Currently Loaded Modules:
  1) gmp/6.3.0-gt7twni (H)   2) gcc/15.1.0

  Where:
   H:  Hidden Module
  ```

- **unload loaded modules**:

  ```
  $ module unload gcc
  ```

- **unload ALL** of the loaded modules:

  ```
  $ module purge
  ```

​

## Private Modules

Users can create [Private Modules](../examples/privatemodule) for the software installed locally.
