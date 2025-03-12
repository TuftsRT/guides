# HPC Storage

The HPC cluster contains a 3PB shared storage system that is available on all the compute nodes.

```{important}
   We will be phasing in charges for storage usage over 10TB. Charges begin July 1, 2026 , for details see the full
   announcement at [RT Announcements](https://it.tufts.edu/research-technology/announcements).
```

```{important}
   No restricted data is allowed on Tufts HPC cluster.
```

The shared storage system is attached to all compute nodes. This means that once you copy your data onto the system, it is available in the same location on any node you connect to.

All users have a home folder at the path `/cluster/home/your_utln` this is where you should store your files. Some users working with a researcher or lab will also have access to a project folder, the path to this will be `/cluster/tufts/your_lab_name/`.

## Storage Usage Limits aka Quotas

```{hint}
   Utilize [hpctools - Tufts HPC Helper Tool](../examples/hpctools) to check storage usage and quota.
```

### Home Directory

Be aware! Your Home Directory (**30GB**, fixed) should be `/cluster/home/your_utln`

If you are not sure how much storage you have used in your home directory, you can find out the information with [hpctools - Tufts HPC Helper Tool](../examples/hpctools) option `6` from command line interface.

### Lab Research Project Storage

Your research project storage (from **50GB and up**) path should be `/cluster/tufts/your_lab_name/`, and each member of the lab group has a dedicated directory `/cluster/tufts/your_lab_name/your_utln`.

See your **research project storage quota** by running the following command from **any node on the cluster**:

```
$ df -h /cluster/tufts/your_lab_name
```

Or utilizing [hpctools - Tufts HPC Helper Tool](../examples/hpctools) option `5` from command line interface.

**NOTE:** Accessing your research project storage space for the **first time** in your current session, please make sure you type out the **FULL PATH** to the top lab directory `/cluster/tufts/your_lab_name/`.

**Need access to a specific lab research project storage on HPC cluster?** Submit a [Cluster Storage Request](https://tufts.qualtrics.com/jfe/form/SV_5bUmpFT0IXeyEfj), the link is also available on [Research Technology website](https://it.tufts.edu/high-performance-computing).

```{gallery-grid}
---
grid-columns: 1
---
- header: "{fas}`book;pst-color-primary` HPC Cluster Storage Request"
  content: "Request research storage space (new and increase) or request access to storage space on HPC cluster."
  link: "../storage-request.html"

- header: "{fas}`book;pst-color-primary` HPC Cluster File Transfers"
  content: "How to transfer files to and from Tufts HPC cluster."
  link: "40-xfer.html"

- header: "{fas}`book;pst-color-primary` File Recovery on Tufts HPC Cluster"
  content: "How to recovery recently deleted files on Tufts HPC cluster using snapshots."
  link: "../examples/snapshots.html"
```
