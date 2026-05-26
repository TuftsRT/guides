# Tier 1 RStore Storage

Tier 1 RStore is high-speed, networked research storage designed for individual researchers or groups who need fast, reliable storage for research instruments, active collaboration and file sharing. Unlike Tier 1 HPC, RStore is not mounted on the HPC compute nodes and is intended for general research use outside of the cluster.

```{important}
Charges of $85/TB/year apply to storage usage exceeding the 10 TB base allocation, effective July 1, 2026.
See [RT Announcements](https://it.tufts.edu/research-technology/announcements) for full details.
```

## Storage Details

| Detail | Value |
|--------|-------|
| Technology | Vast, SSD |
| Starting quota | 100 GB minimum |
| Base allocation | 10 TB at no cost (if eligible faculty member)|
| Cost beyond base | $85/TB/year (effective July 1, 2026) |

## Prerequisites

- Active Tufts account
- Faculty sponsorship for research storage share
- Off-campus access requires [Tufts VPN](https://access.tufts.edu/vpn)

## Requesting Storage

Submit a [Research Storage Request](https://tufts.qualtrics.com/jfe/form/SV_5bUmpFT0IXeyEfj) to:

- Request a new RStore share
- Increase an existing quota

Only Tufts University faculty can be the PI of a research storage share. Use the [Data Storage Finder](https://access.tufts.edu/data-finder) to confirm RStore is the right tier for your needs.


## Access Methods

### Network Drive (SMB/NFS)

RStore is typically connected directly to end user workstations or laptops. It can be mounted as a network drive on Windows or Mac for direct file access over the Tufts network. Off-campus access requires [Tufts VPN](https://access.tufts.edu/vpn).

#### Mac

1. Open **Finder** and press **Command+K**.
2. Enter the path to your RStore share:
   ```
   smb://rstore.it.tufts.edu/RStoreDriveName
   ```
3. Click **Connect**.

#### Windows

1. Open **File Explorer** and right-click **This PC**.
2. Select **Map network drive...**.
3. In the **Drive** dropdown, select a drive letter (R, S, T, etc.).
4. In the **Folder** text box, enter the path to your RStore share:
   ```
   \\rstore.it.tufts.edu\RStoreDriveName
   ```
5. For **Tufts-owned computers**, click **Finish**.

   For **personal computers**, check **Connect using different credentials** before clicking Finish, then enter your Tufts credentials at the prompt:
   - **Username:** `TUFTS\Tufts_Username`
   - **Password:** `Tufts_Password`

   Confirm the domain shown reads **TUFTS**, then click **OK**.



