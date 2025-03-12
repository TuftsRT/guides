# Private Modules

Users on Tufts HPC cluster can create their own module files and load their own private modules to make using locally installed software easier.

## 1. Install your software

After installation, make sure you know the path to the software. Especially the`bin` `lib` or `lib64` `include` and `man` if there are any.

## 2. Create your module file

The **location** of the module file is very important! If it's misplaced, your module won't show up.

Please put your own module file in your home directory, under `privatemodules`.

Each software should have its own directory, such as `mysoftware`.

And for each version of the software, there should be a module file, such as `v1.0.0`.

The overall structure should be:

`/cluster/home/$USER/privatemodules/mysoftware/v1.0.0`

Now let's work on the content of the file.

You can reference/borrow existing module files on the cluster.

For simple software, you can use "cmake/3.18" as a template.

`$ module av cmake`

```
------------------------------------------------------------------------------------------------------ /opt/shared/Modules/modulefiles-rhel6 ------------------------------------------------------------------------------------------------------
cmake/2.8      cmake/2.8.11.2 cmake/3.2.1    cmake/3.4.3

--------------------------------------------------------------------------------------------------------- /cluster/tufts/hpc/tools/module ---------------------------------------------------------------------------------------------------------
cmake/3.18     cmake/3.23_gui
```

To take a look at the module file:

`$ cat /cluster/tufts/hpc/tools/module/cmake/3.18`

```
#%Module -*- tcl -*-
#
# cmake 3.18
#
#

set ver 3.18

proc ModulesHelp { } {
  global ver
  puts stderr "\tThis module sets up the cmake 3.18 environment installed by DM in hpc_tools, version $ver "
}

module-whatis   "load cmake 3.18"


prepend-path PATH /cluster/tufts/hpc/tools/spack/linux-rhel7-ivybridge/gcc-8.4.0/cmake-3.18.2-r2hx6rewcozqvda5gfo7rouhzlqy2k6t/bin

#
# appended log section
#

if {[module-info mode "load"]} {
  global env
  if {[info exists env(USER)]} {
    set the_user [lindex [array get env USER] 1]
  } else {
    set the_user "foo"
  }
  system [concat "logger environment-modules" [module-info name] $the_user ]
}

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
------------- /cluster/home/$USER/privatemodules ---------------------------------------------------------------------------------------------------------
mysoftware/v1.0.0
```

You can load it just like other modules on the cluster:

`$ module load my software/v1.0.0`

:)
