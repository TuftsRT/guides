# Tier 1 HPC Storage

Tier 1 HPC is a multi-petabyte, all-flash storage system attached to the Tufts HPC cluster. It is available on all compute and login nodes, meaning data stored here is accessible from any node without copying.

```{important}
Charges of $85/TB/year apply to storage usage exceeding the 10 TB base allocation, effective July 1, 2026.
See [RT Announcements](https://it.tufts.edu/research-technology/announcements) for full details.
```

## Storage Locations

| Location            | Path                            | Default Quota |
| ------------------- | ------------------------------- | ------------- |
| Home directory      | `/cluster/home/your_utln`       | 30 GB (fixed) |
| Lab/project storage | `/cluster/tufts/your_lab_name/` | 50 GB and up  |

## Storage Details

| Detail           | Value                                         |
| ---------------- | --------------------------------------------- |
| Technology       | Vast, SSD, Parallel Filesystem                |
| Starting quota   | 100 GB minimum                                |
| Base allocation  | 10 TB at no cost (if eligible faculty member) |
| Cost beyond base | \$85/TB/year (effective July 1, 2026)         |

## Prerequisites

- Active Tufts HPC cluster account — see [HPC Account Request](../../hpc/access/10-account-request.md)
- Off-campus access requires [Tufts VPN](https://access.tufts.edu/vpn)

## Requesting Storage

- **New storage share:** Submit a [Cluster Storage Request](https://tufts.qualtrics.com/jfe/form/SV_5bUmpFT0IXeyEfj)
- **Increase existing quota:** Submit the same form; increases over 1 TB may require a consultation
- No more than 5 TB increase per 12-month interval per current policy

## Access Methods

### OnDemand (Web Browser)

The simplest way to access and manage files without any software installation.

1. Go to [Tufts HPC OnDemand](https://ondemand-prod.pax.tufts.edu/)
1. Log in with your Tufts credentials
1. Click **Files** in the top menu to browse your home directory
1. Use **Change directory** to navigate to your lab folder — type the full path (e.g., `/cluster/tufts/your_lab_name/`)
1. Use **Upload** / **Download** buttons to transfer files

```{note}
OnDemand file transfers are limited to files **under 976 MB**. Use SCP, rsync, or Globus for larger files.
```

### SSH / Command Line

Connect to the cluster login node and access storage directly:

```
ssh your_utln@login-prod.pax.tufts.edu
```

Your home directory is available immediately. Navigate to your lab storage:

```
cd /cluster/tufts/your_lab_name/
```

### SCP / rsync (File Transfer)

Run these commands from your **local machine**. Replace `your_utln` and paths accordingly.

**Upload to cluster:**

```
scp local_file your_utln@login-prod.pax.tufts.edu:/cluster/home/your_utln/
rsync -azP local_dir/ your_utln@login-prod.pax.tufts.edu:/cluster/tufts/your_lab_name/
```

**Download from cluster:**

```
scp your_utln@login-prod.pax.tufts.edu:/cluster/home/your_utln/file ./
rsync -azP your_utln@login-prod.pax.tufts.edu:/cluster/tufts/your_lab_name/ local_dir/
```

Use `rsync` for large or long-running transfers — it can resume if interrupted.

Recommended GUI clients: [SecureCRT](https://access.tufts.edu/securecrt), [FileZilla](https://filezilla-project.org/), [Cyberduck](https://cyberduck.io/), [WinSCP](https://winscp.net/) (Windows)

### Globus

Globus is recommended for large or automated transfers. Tufts provides a site collection for HPC cluster storage.

- Collection name: **Tufts HPC Cluster Storage**
- Provides access to `/cluster/home/your_utln` and `/cluster/tufts/`

See the [Globus guide](../globus/index.md) for setup and transfer instructions.

## Checking Storage Usage

Use `hpctools` on any cluster node via SSH Shell. Go to [Tufts HPC Cluster OnDemand](https://ondemand-prod.pax.tufts.edu/) click , **Clusters** then **Tufts HPC Shell Access** or log in to the cluster through command line with SSH.

```
module load hpctools
hpctools
```

- Option **5** — check lab/project storage quota
- Option **7** — check home directory usage

Or check project quota directly:

```
df -H /cluster/tufts/your_lab_name
```
