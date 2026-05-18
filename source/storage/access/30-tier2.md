# Tier 2 Research Storage

Tier 2 is a reliable, scalable storage environment for inactive research data that no longer needs to reside on high-performance Tier 1 systems. It is connected to Tier 1 (HPC and RStore) via high speed (100Gbps) network links, making data migration straightforward and fast.

```{important}
Charges of $35/TB/year apply to storage usage exceeding the 10 TB base allocation, effective July 1, 2026.
See [RT Announcements](https://it.tufts.edu/research-technology/announcements) for full details.
```

## Storage Details

| Detail | Value |
|--------|-------|
| Path | `/cluster/tier2/projectname` |
| Starting quota | 100 GB minimum |
| Base allocation | 10 TB at no cost (if eligible faculty member)|
| Cost beyond base | $35/TB/year (effective July 1, 2026) |

## Prerequisites

- Active Tufts account and research storage allocation
- A [Globus account](../globus/42-globus-account-setup.md) linked to your Tufts credentials

## Access Method — Globus

Tier 2 storage can be accessed through **Globus**. Tufts Research Technology provides a site collection for Tier 2.

1. Log in to [Globus](https://www.globus.org/) with your Tufts credentials
2. In the **File Manager**, click the **Collection** field and search for **"Tufts"**
3. Select the **Tufts Tier 2 Storage** collection
4. Authenticate with your Tufts credentials when prompted
5. Your `/cluster/tier2/projectname` directory will be available to browse and transfer

See the full [Globus guide](../globus/index.md) for setup, transferring files between endpoints, and monitoring transfers.

## Moving Data Between Tiers

A common workflow is to migrate completed or infrequently accessed data from Tier 1 HPC or RStore to Tier 2:

### To Transfer using Globus
1. Set up both a Tufts Tier 1 endpoint and a Tufts Tier 2 endpoint in Globus File Manager
2. Navigate to the source folder on the Tier 1 collection
3. Navigate to the destination folder on the Tier 2 collection
4. Select files and click **Transfer or Sync to**

### To Transfer via Tufts HPC Cluster login node

Tier 2 is mounted on the HPC login nodes, so you can move data directly from the command line without Globus.

1. SSH into the login node:

```
ssh your_utln@login-prod.pax.tufts.edu
```

2. Move data from Tier 1 to Tier 2:

```
mv /cluster/tufts/projectname /cluster/tier2/projectname
```

To move a specific subdirectory or file:

```
mv /cluster/tufts/projectname/subdir /cluster/tier2/projectname/subdir
```

```{note}
`mv` is preferred over `cp` when migrating data to Tier 2 to avoid duplicate storage charges. If you need to keep a copy on Tier 1, use `cp -r` instead.
```



## Requesting Storage

Submit a [Research Storage Request](https://tufts.qualtrics.com/jfe/form/SV_5bUmpFT0IXeyEfj) to request a new Tier 2 allocation or increase an existing one. Use the [Data Storage Finder](https://access.tufts.edu/data-finder) to confirm Tier 2 is appropriate for your data.
