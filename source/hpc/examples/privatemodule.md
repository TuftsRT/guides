# Private Modules

Users on Tufts HPC cluster can create their own module files and load their own private modules to make using locally installed software easier.

## 1. Install your software

After installation, make sure you know the path to the software. Especially the`bin` `lib` or `lib64` `include` and `man` if there are any.

## 2. Create your module file

The **location** of the module file is very important! If it's misplaced, your module won't show up.

Please put your own module file in your home directory, under `privatemodules`.

Each software should have its own directory, such as `mysoftware`.

And for each version of the software, there should be a module file, such as `1.0.0.lua`.

The overall structure should be:

`/cluster/home/$USER/privatemodules/mysoftware/1.0.0.lua`

Now let's work on the content of the file.

You can refer to the following `cmake` modulefile located at `/cluster/tufts/apps/modules/9/x86_64/gcc/15.1.0/cmake/3.31.6.lua`.

```
-- -*- lua -*---
-- cmake@3.31.6~doc+ncurses+ownlibs+qtgui build_system=generic build_type=Release arch=linux-rocky9-broadwell/htsgau3
--

whatis([[Name : cmake]])
whatis([[Version : 3.31.6]])
whatis([[Short description : A cross-platform, open-source build system. CMake is a family of tools designed to build, test and package software. ]])

help([[A cross-platform, open-source build system. CMake is a family of tools
designed to build, test and package software.]])


depends_on("gmake/4.4.1-yi34qio")
depends_on("qt/5.15.12-qkjcg5f")

local modroot="/cluster/tufts/apps/spack/9/x86_64/apps/linux-broadwell/cmake-3.31.6-htsgau3kr7aenjwhvsdbkadkwjrfgulh"
prepend_path("PATH", modroot.."/bin", ":")
prepend_path("ACLOCAL_PATH", modroot.."/share/aclocal", ":")
prepend_path("CMAKE_PREFIX_PATH", modroot.."/.", ":")
```

Most of the content is informational.

The most important function part of the module file would be the `prepend-path` line, it adds the cmake `bin` directory path to the environment variable `PATH` so that the system will be able to find the cmake executables.

If your software requires other environment variables such as `LD_LIBRARY_PATH` or `CPATH` to be set to function, then you should `prepend-path` those variables as well with the appropriate directories.

Once your module file is ready. Save it.

## 3. Load Your Own Module

In order to enable the usage of private module, `use.own` module needs to be loaded first:

` $ module load use.own`

Now you can check your own module.

`$ module av mysoftware`

And it should show up as:

```
------------- /cluster/home/$USER/privatemodules ------------------
mysoftware/1.0.0
```

You can load it just like other modules on the cluster:

`$ module load my software/1.0.0`

:)
