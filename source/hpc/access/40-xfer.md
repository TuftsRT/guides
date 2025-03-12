# File Transfers

```{important}
   [Globus](Globus/index.md) is highly recommended for File Transfers on Tufts HPC Cluster.
```

```{important}
   VPN - Off-campus access from non-Tufts Network please connect to [Tufts VPN](https://access.tufts.edu/vpn).
```

## Prerequisites

- Active Cluster Account
- VPN
- 2FA

There are several methods available to copy your data to and from the HPC cluster. These include SCP/SFTP, OnDemand and Globus.

```{important}
   File Transfer Hostname xfer.pax.tufts.edu (formerly xfer.cluster.tufts.edu)
```

- File Transfer Protocol

  - SCP(Secure Copy Protocol)
  - SFTP(SSH File Transfer Protocol)
  - rsync over SSH

- Recommended software

  > - Cross-Platform: **[SecureCRT](https://access.tufts.edu/securecrt)**
  > - Windows Only: **[WinSCP](https://winscp.net/eng/index.php)**
  > - Cross-Platform: **[FileZilla](https://filezilla-project.org/)**
  > - Cross-Platform: **[Cyberduck](https://cyberduck.io/)**
  > - Native: Windows 11, Mac and Linux all include a command line SSH clients in the OS. This `scp` command can be used for quick connections without installing additional software.
  > - Web Interface:**[OnDemand](https://ondemand.pax.tufts.edu)** (Only for transferring files size **less than 976MB per file.**)

## OnDemand

- Go to OnDemand:

  **[https://ondemand.pax.tufts.edu/](https://ondemand.pax.tufts.edu/)**

- Under **`Files`**, using the **`Upload`** or **`Download`** buttons to transfer. Make sure you navigate to the destination/source directory on cluster using the **`Go To`** button before transferring files.

  <img src="https://raw.githubusercontent.com/DelilahYM/ImageHost/master/ondemand_homedir.png" alt="Core-Node" width=70%>

  You may not see your lab folder listed under `/cluster/tufts`, **DO** use `Go To` button and **TYPE** out the **FULL** path.

## Command Line - scp & rsync

> Hostname\*\* for file transfer: xfer.cluster.tufts.edu
>
> NOTE:
>
> - **Local_Path** : Path to your files or directory on your local computer
> - **Cluster_Path** :Path to your files or directory on the cluster
>   - Home Directory: _/cluster/home/your_utln/your_folder_
>   - Research Project Storage Space Directory: _/cluster/tufts/your_lab_name/your_utln/your_folder_

\*\*_Execute from your local machine terminal._ \*\*

General Format:

`$ scp From_Path To_Path`

`$ rsync From_Path To_Pat`

**_NOTE: If you are transferring very large files that could take hours to finish, we would suggest using `rsync` as it has ability to restart from where it left if interrupted._**

Example:

**File** Transfer with `scp`or `rsync`:

- Download from cluster

`$ scp your_utln@xfer.cluster.tufts.edu:Cluster_Path Local_Path  `

`$ rsync your_utln@xfer.cluster.tufts.edu:Cluster_Path Local_Path`

- Upload to cluster

`$ scp Local_Path your_utln@xfer.cluster.tufts.edu:Cluster_Path`

`$ rsync Local_Path your_utln@xfer.cluster.tufts.edu:Cluster_Path`

**Directory** Transfer with `scp` or `rsync`:

- Download from cluster

`$ scp -r your_utln@xfer.cluster.tufts.edu:Cluster_Path Local_Path  `

`$ rsync -azP your_utln@xfer.cluster.tufts.edu:Cluster_Path Local_Path`

- Upload to cluster

`$ scp -r Local_Path your_utln@xfer.cluster.tufts.edu:Cluster_Path`

`$ rsync -azP Local_Path your_utln@xfer.cluster.tufts.edu:Cluster_Path`
