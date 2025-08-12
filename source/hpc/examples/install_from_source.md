# Building Software from Source Code on the Cluster

Users can also compile/install oftware from source code on Tufts HPC. This guide will walk you through the essential tools and steps you'll need to successfully build and install applications.

## ⚙️ Compilers: The Core of the Process
A compiler is a critical program that translates your application's source code (e.g., C++ or Fortran) into machine-readable instructions. The choice of compiler can significantly impact your application's performance.

On our cluster, you have access to two powerful compiler suites:
- GNU Compiler Collection (GCC): A versatile and widely used open-source compiler suite. It's an excellent general-purpose option and the default for many applications.

- Intel oneAPI Compilers: A suite of compilers specifically optimized for Intel processors, which are used on our cluster. Using this suite can lead to substantial performance gains for scientific and HPC applications.
You can check available compiler versions using the module avail command:

### Check available GCC versions
```
$ module avail gcc
gcc/12.4.0    gcc/15.1.0 (D)
```

###  Check available Intel oneAPI compiler versions
```
$ module avail intel-oneapi-compilers
intel-oneapi-compilers/2025.1.1
```

## 🔨 GNU Make: Automating the Build Process
GNU Make is a program that automates the process of building executables from source code. It uses a file called Makefile (or makefile) which contains a set of rules for compiling the program.

### General Steps for Installation with Make
- Unpack the source code: Download and extract the compressed source files.
- Load the compiler: Use module load to activate the correct compiler (e.g., GCC).
- Configure (if applicable): Some packages have a configure script that checks your system and creates a custom Makefile.
- Build the program: Run the make command to compile the source code.
- Install the software (optional): Run make install to move the compiled files to a designated location. **Do not use sudo make install, as you don't have administrative privileges on the cluster**.

`Tip`: By default, make install tries to install into system directories like `/usr/local` where you lack permissions. To install to a custom directory (like your home directory), you must configure the installation path. For packages with a configure script, use the `--prefi`x option: `./configure --prefix=$HOME/my-app`.

### Example: Installing BWA with Make
BWA is a popular tool for mapping DNA sequences. This example shows a typical workflow for installing software with make.
Set up your installation directory: It's a good practice to create a dedicated directory for your custom software.
```
$ cd $HOME
$ mkdir apps
$ cd apps
```

Load the required compiler: BWA requires a C compiler, so we'll load the latest GCC module.

```
$ module load gcc/15.1.0
```
Download the source code: Clone the BWA repository from GitHub.
```
$ git clone https://github.com/lh3/bwa.git
```
Build the application: Navigate into the new directory and run make.
```
$ cd bwa
$ make
```

Add the new program to your PATH: The compiled bwa executable is now in the bwa directory. To run it from any location, you need to add this directory to your PATH environment variable.
```
$ export PATH=$PATH:$HOME/apps/bwa
```
Verify the installation: You can now run the bwa command to see if it works.
```
$ bwa
Program: bwa (alignment via Burrows-Wheeler transformation)
```

## 🧱 CMake: A Higher-Level Build System
CMake is a popular cross-platform build system generator. Instead of compiling software directly, it generates native build scripts (like Makefiles) that make can then use.

### Example: Installing regtools with CMake
Some applications, like regtools, require cmake to generate the Makefile before make can be used.
Load the necessary modules: Always load both the compiler (e.g., gcc) and the latest cmake version.
```
$ module load gcc/15.1.0
$ module load cmake/3.31.6
```

Clone the source code:
```
$ cd $HOME/apps
$ git clone https://github.com/griffithlab/regtools
```

Build the application: It's a best practice to create a separate build directory to keep your source code clean.
```
$ cd regtools/
$ mkdir build
$ cd build/
```
Run cmake to generate the Makefile, then use make
```
$ cmake ..
$ make
```

### Install to a custom prefix

For applications that use `make install`, you must specify the installation directory with the `-DCMAKE_INSTALL_PREFIX` option during the cmake step.
This is a common pattern for installing CMake-based applications.
```
$ cd $HOME/apps/regtools
$ mkdir build
$ cd build

$ cmake -DCMAKE_INSTALL_PREFIX=$HOME/apps/regtools-install ..
$ make
$ make install
```

This command will build the program and then move the final files into a new directory named regtools-install. You can then add the bin directory from this new path to your `$PATH` variable.
