---
tags: hpc
---

# 2025 Cluster Upgrade

```{attention}
The upgraded cluster is currently only accessible to approved early adopters.  
You may request to join our early adopter program at https://tufts.qualtrics.com/jfe/form/SV_08IS0n1YSTR6KRU .
```

In the summer of 2025 TTS Research Technology launched a new version of the pax cluster with a newer OS version.  This system shares a 
file system with the existing cluster so all your files and data are accessible in both versions.

## Migration Plan
New users will be encouraged to start using this version when they receive their accounts.  Existing HPC users are 
encouraged to try this version as time allows.  As the average compute load on the upgraded cluster increases we will 
move public partition nodes from the current cluster to this version.

TTS Research Technology will work with contribute node owners to confirm their software works on the upgraded cluster, 
and determine a timeframe to move their nodes as the labs convenience.

## Upgraded cluster access
To access the upgraded cluster

* SSH: login-prod.pax.tufts.edu
* Web: https://ondemand-prod.pax.tufts.edu/

## Notable differences

* Operating systems is Linux Rocky 9 (RHEL9 compatible) providing the latest system libraries including; glibc 
* Upgraded cluster is running Nvidia Driver Version: 575.57.08 with CUDA Version: 12.9

## Modules
The upgraded cluster continues to use the modules system to load software.  However the software available is 
probably a newer version. All software names have been standardized in lower case.

* The old RHEL7 cluster software is available if you need an older version of software. Run the command `module load 
  modtree/deprecated` to switch to the old modules.  These modules may or may not work on the new cluster, and you should 
  work to migrate to more recent versions of your software.  To switch back to the current modules, run `module unload modtree/deprecated`.  If you need 
  assistance with this contact RT.

## Known Issues
The following is a list of known issues that are currently being worked on.

* Broken open ondemand applications
  * Matlab server
  * Fiji, and any other VNC based applications

## Upgraded cluster resources

New GPU resources are available on this cluster. 
* H200 GPUs , request using `--gres=gpu:h200:1` 


