## 2025 Cluster Upgrade

```{attention}
The upgraded cluster is currently only accessible to approved early adopters.
You may request to join our [early adopter program](https://tufts.qualtrics.com/jfe/form/SV_08IS0n1YSTR6KRU).
```
In the summer of 2025, TTS Research Technology upgraded the pax cluster to a new version with an updated operating system. The new cluster shares a file system with the existing one, so all your files and data remain accessible on both systems.

![](img/os.png){width=700}

##### Why do we need to upgrade cluster's operation system?
- **End of Life**  
  RHEL7 reached end-of-life in June 2024 and no longer receives official security patches or bug fixes. RHEL8 will also phase out sooner than Rocky 9.6’s supported lifecycle.  

- **Outdated compilers and libraries**  
  Many modern HPC applications cannot be built or run on the old compiler toolchains and libraries available in RHEL7.  

- **Incompatibility with new hardware**  
  Older kernels from RHEL7/8 lack support for modern GPUs, CPUs, interconnects, and storage technologies used in today’s HPC clusters.  

- **Security and compliance risks**  
  Running an unsupported OS leaves the system exposed to known vulnerabilities and makes it difficult to comply with security requirements.  

- **Vendor and community support has moved on**  
  Most open-source and commercial software providers no longer test or release updates for RHEL7, meaning bugs or failures will not be addressed.  


##### Benefits of the new Rocky Linux 9.6 operation system
- **Long-term support**  
  Rocky Linux 9.6 is a modern, community-supported enterprise OS with long lifecycle support, ensuring stability and security updates for years to come.  

- **Improved security**  
  The updated OS includes the latest security patches, kernel improvements, and system hardening features to better protect user data and workloads.  

- **Modern software ecosystem**  
  Many scientific applications and libraries no longer support RHEL7/8. Rocky Linux 9.6 provides access to newer compilers, MPI stacks, and optimized math libraries that improve performance and compatibility.  

- **Better hardware support**  
  The newer kernel and drivers provide improved support for GPUs, high-speed interconnects, and other modern HPC hardware.  

- **Performance improvements**  
  Applications rebuilt on Rocky Linux 9.6 often see better runtime performance due to compiler optimizations and updated toolchains.  


##### Impacts to Users

The move to Rocky Linux 9.6 will bring noticeable changes to the software environment. Users should be aware of the following:

- **Recompiled software stack**  
  All applications have been rebuilt on the new OS. You may see differences in versions and behavior or performance compared to the old RHEL7/8 versions.  

- **Legacy modules**  
  Older RHEL6/7/8 software stacks are still available using:  
  - `module load modtree/deprecated`  
  However, these modules may not function reliably on Rocky 9.6 and will not be maintained long-term.  

- **Testing Your Job Scripts is Recommended**  
  You'll need to test any applications used in your current job scripts on the new clusters, as they may not work without an update.
  - **Update Software Modules**: If newer software modules are available, we recommend updating your job scripts to use them.
  - **Rebuild Self-Built Software**: If you're using self-built software, you may need to rebuild it on the new cluster.
 If you have any issues rebuilding your software or migrating your workflows, please contact us at, please contact us at [tts-research@tufts.edu](mailto:tts-research@tufts.edu).

- **Migration required**  
  While the legacy stacks are available temporarily, users should plan to migrate to the Rocky 9.6 modules as soon as possible. Future support will focus exclusively on the new environment.  

#### Upgraded cluster access

To access the upgraded cluster

- SSH: ssh yourUTLN@login-prod.pax.tufts.edu
- New Open OnDemand: https://ondemand-prod.pax.tufts.edu/

#### Migration Timeline
- **New Users**: New users will be encouraged to start using the new cluster when they receive their accounts.

- **Existing Users**: Existing HPC users are encouraged to try the new cluster as time allows. As the average compute load on the upgraded cluster increases, we will move public partition nodes from the current cluster to this new cluster.

- **Contribute Node Owners**: TTS Research Technology will work with contribute node owners to confirm their software works on the upgraded cluster and determine a timeframe to move their nodes at the lab's convenience.


#### New hardwares

**NVIDIA H200 GPU**

[NVIDIA H200](https://www.nvidia.com/en-us/data-center/h200/) is a powerhouse for large-scale AI and HPC. It provides a massive boost in memory capacity and bandwidth, which is crucial for training and running large language models and complex simulations.

  - Compute Capability: 9.0 (Hopper Architecture)
  - GPU Memory: 141 GB HBM3e
  - Requesting a GPU: To submit a job to a node with an H200 GPU, use the following SLURM directive in your job script:
    ```
    #SBATCH --partition=gpu
    #SBATCH --gres=gpu:h200:1
    ```

**NVIDIA L40s GPU (coming soon)**

[NVIDIA L40s]((https://www.nvidia.com/en-us/data-center/l40s/)) is a versatile, universal GPU designed for a wide array of data center workloads. It excels at multi-modal generative AI, 3D graphics, and video processing, offering a balance of compute and visualization capabilities.

  - Compute Capability: 8.9 (Ada Lovelace Architecture)
  - GPU Memory: 48 GB GDDR6
  - Requesting a GPU: To submit a job to a node with an L40s GPU, use the following SLURM directive in your job script:
    ```
    #SBATCH --partition=gpu
    #SBATCH --gres=gpu:l40s:1
    ```

#### Applications and servers

```{gallery-grid}
---
grid-columns: 1
---
- header: "{fas}`flask` Applications"
  content: "A wide selection of pre-built scientific applications that users can easily access and load as modules, eliminating the need for self-installation."
  link: "module.html"

- header: "{fas}`laptop-code` Open OnDemand"
  content: "A new Open OnDemand server that provides a more stable, reliable platform and a modern, intuitive interface with features designed to improve your workflow."
  link: "ood.html"
```
