# File Transfers

```{important}
   [Globus](globus/index) is highly recommended for File Transfers on Tufts HPC Cluster.
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
   File Transfer hosts (xfer) are being retired in favor of using the login nodes directly for SCP/SFTP or rsync.
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
  > - Web Interface:**[OnDemand](https://ondemand-prod.pax.tufts.edu)** (Only for transferring files size **less than 976MB per file.**)

## OnDemand

- Go to OnDemand:

  **[https://ondemand-prod.pax.tufts.edu/](https://ondemand-prod.pax.tufts.edu/)**

- **`Home Directory`** Under **`Files`**, using the **`Upload`** or **`Download`** buttons to transfer. Make sure you navigate to the destination/source directory on cluster using the **`Change directory`** button before transferring files.

![ood-homedir](../assets/el9/newondemand-homedir.png)

  You may not see your lab folder listed under `/cluster/tufts`, **DO** use `Change directory` button and **TYPE** out the **FULL** path.

## Command Line - scp & rsync

> Hostname for file transfer: **login-prod.pax.tufts.edu**
>
> NOTE:
>
> - **Local_Path** : Path to your files or directory on your local computer
> - **Cluster_Path** :Path to your files or directory on the cluster
>   - Home Directory: _/cluster/home/your_utln/your_folder_
>   - Research Project Storage Space Directory: _/cluster/tufts/your_lab_name/your_utln/your_folder_

The following commands are execute from your local machine terminal.

General Format:

`$ scp From_Path To_Path`

`$ rsync From_Path To_Pat`

**_NOTE: If you are transferring very large files that could take hours to finish, we would suggest using `rsync` as it has ability to restart from where it left if interrupted._**

Example:

**File** Transfer with `scp`or `rsync`:

- Download from cluster

`$ scp your_utln@login-prod.pax.tufts.edu:Cluster_Path Local_Path  `

`$ rsync your_utln@login-prod.pax.tufts.edu:Cluster_Path Local_Path`

- Upload to cluster

`$ scp Local_Path your_utln@login-prod.pax.tufts.edu:Cluster_Path`

`$ rsync Local_Path your_utln@login-prod.pax.tufts.edu:Cluster_Path`

**Directory** Transfer with `scp` or `rsync`:

- Download from cluster

`$ scp -r your_utln@login-prod.pax.tufts.edu:Cluster_Path Local_Path  `

`$ rsync -azP your_utln@login-prod.pax.tufts.edu:Cluster_Path Local_Path`

- Upload to cluster

`$ scp -r Local_Path your_utln@login-prod.pax.tufts.edu:Cluster_Path`

`$ rsync -azP Local_Path your_utln@login-prod.pax.tufts.edu:Cluster_Path`
