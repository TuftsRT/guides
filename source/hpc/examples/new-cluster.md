---
tags: hpc
---

# HPC Cluster Upgrade

```{attention}
The upgraded cluster is currently only accessible to approved early adopters.
You may request to join our early adopter program at https://tufts.qualtrics.com/jfe/form/SV_08IS0n1YSTR6KRU .
```

In the summer of 2025 TTS Research Technology launched a new version of the pax cluster with a newer OS version. This system shares a
file system with the existing cluster so all your files and data are accessible in both versions.

- Newer OS provides longer product life and support
- Strengthens security
- Up-to-date drivers, compilers and libraries provide a reliable foundation for software installation and execution
- Updated system unlocks the power of next-generation GPUs, CPUs, interconnects, and storage with seamless compatibility
- Establishing a strong foundation for future enhancements

## Migration Plan

New users will be encouraged to start using this version when they receive their accounts. Existing HPC users are
encouraged to try this version as time allows. As the average compute load on the upgraded cluster increases we will
move public partition nodes from the current cluster to this version.

TTS Research Technology will work with contribute node owners to confirm their software works on the upgraded cluster,
and determine a timeframe to move their nodes as the labs convenience.

## Compute nodes

Public partition nodes are being moved to the new cluster throughout the spring of 2026. Most public resources will be moved by the end of May. Please test your software and workflows, and move your usage to the upgraded cluster ASAP. Significant planned migration dates are listed below.

- 32 A100 GPUs migrating on 3/1/2026
- All public GPUs migrated by 5/31/2026

TTS Research Technology will work with contribute node owners to confirm their labs software works on the upgraded cluster, and determine a timeframe to move their nodes at the labs convenience.

## Software / Service retirement

As part of this upgrade several legacy software and services are being retired. We have identified a modern replacement for each. If you have concerns or questions about the replacements, or timeline please contact us.

| Software being retired | Replacement    | Retirement Date |
| ---------------------- | -------------- | --------------- |
| Galaxy                 | NF-Core        | 2/1/2026        |
| FastX                  | HPC Desktop    | 4/1/2026        |
| Flamenco               | Blender Render | 6/30/2026       |

## Upgraded cluster access

To access the upgraded cluster

- SSH: login-prod.pax.tufts.edu
- Web: https://ondemand-prod.pax.tufts.edu/

## Notable differences

- Operating systems is Rocky Linux 9 (RHEL9 compatible) providing the latest system libraries including: glibc
- Hyper-Threading has been disabled on all compute nodes in the upgraded cluster
- Upgraded cluster is running Nvidia Driver Version: 575.57.08 with CUDA Version: 12.9
- No GPUs in the preempt partition, use the regular gpu partition if you need GPUs
- MPI partition has been retired, use the batch partition

## Modules

The upgraded cluster continues to use the modules system to load software. However the software available is
probably a newer version. All software names have been standardized in lower case.

- The old RHEL7 cluster software is available if you need an older version of software. Run the command `module load  modtree/deprecated` to switch to the old modules. These modules may or may not work on the new cluster, and you should
  work to migrate to more recent versions of your software. To switch back to the current modules, run `module unload modtree/deprecated`. If you need
  assistance with this contact RT.

## Known Issues

The following is a list of known issues that are currently being worked on.

- Broken open ondemand applications
  - Matlab server
  - Fiji, and any other VNC based applications

## Upgraded cluster resources

New GPU resources are available on this cluster.

- H200 GPUs , request using `--gres=gpu:h200:1`
- L40s GPUs , request using `--gres=gpu:l40s`
