---
tags: hpc
---

# 2025 Cluster Upgrade

```{attention}
The upgraded cluster is currently only accessible to approved early adopters.  
You may request to join our early adopter program at https://tufts.qualtrics.com/jfe/form/SV_08IS0n1YSTR6KRU .
```

In the summer of 2025 RT launched a new version of the pax cluster with a newer OS version.

As the average compute load on the upgraded cluster increases we will move public partition nodes from the current
cluster.

RT will work with contribute node owners to confirm their software works on the upgraded cluster, and determine
a timeframe that works for the lab to move their nodes.

## Upgraded cluster access
To access the upgraded cluster

* SSH: login-prod.pax.tufts.edu
* Web: https://ondemand-prod.pax.tufts.edu/

## Notable differences

* Upgraded cluster is running Nvidia Driver Version: 575.57.08 with CUDA Version: 12.9

## Modules
The upgraded cluster uses the modules system to load software.  The software available may be newer version, and all 
software names have been standardized in lower case.

* The old RHEL7 cluster software is available if you need an older version of software. Run the command `module load 
  deprecated` to switch to the old modules.  These modules may or may not work on the new cluster, and you should 
  work to migrate to more recent versions of your software.  If you need assistance with this contact RT.


## Upgraded cluster resources

New GPU resources are available on this cluster. 
* H200 GPUs , request using `--gres=gpu:h200:1` 


