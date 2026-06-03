# Tier 3 Research Storage

Tier 3 is a low-cost tape-based storage environment for static, rarely accessed research data. It is operated through the **NESE (Northeast Storage Exchange) Tape** service, maintained by Harvard at MGHPCC. Tier 3 is connected to Tier 2 via high-speed network links for straightforward data migration.

```{important}
Tier 3 has **no free base allocation** — all storage is charged at $12/TB/year, effective July 1, 2026.
See [RT Announcements](https://it.tufts.edu/research-technology/announcements) for full details.
```

```{important}
Tier 3 is available for **HPC storage only** and is intended for long-term retention of data that does not require regular access.
```

## Storage Details

| Detail             | Value                                 |
| ------------------ | ------------------------------------- |
| Technology         | Tape (NESE, hosted at MGHPCC)         |
| Minimum allocation | 5 TB (in 5 TB increments)             |
| Base allocation    | None — all storage is charged         |
| Cost               | \$12/TB/year (effective July 1, 2026) |

## Prerequisites

- Active Tufts account and a Tier 2 or Tier 1 HPC storage allocation
- A [Globus account](../globus/42-globus-account-setup.md) linked to your Tufts credentials

## Requesting Storage

Submit a [Research Storage Request](https://tufts.qualtrics.com/jfe/form/SV_5bUmpFT0IXeyEfj) to request a Tier 3 allocation. Allocations are provisioned in 5 TB increments. Use the [Data Storage Finder](https://access.tufts.edu/data-finder) to confirm Tier 3 is appropriate for your data.

For questions, contact [tts-research@tufts.edu](mailto:tts-research@tufts.edu).

## Access Method — Globus

Tier 3 tape storage is accessed exclusively through **Globus**.

1. Log in to [Globus](https://www.globus.org/) with your Tufts credentials
1. In the **File Manager**, click the **Collection** field and search for **"Tufts"**
1. Select the **Tufts Tier 3 Storage** collection
1. Authenticate with your Tufts credentials when prompted
1. Transfer files to or from your Tier 3 allocation

See the full [Globus guide](../globus/index.md) for setup, transferring files, and monitoring transfers.

```{note}
Tape storage has higher latency than disk-based tiers. Retrieval of data may take longer than transfers between Tier 1 and Tier 2.
```

## Moving Data to Tier 3

Tier 3 is best suited for data that has already been moved to Tier 2 and is no longer accessed regularly. Archiving data to Tier 3 is not a self service operation. If you have data in Tier 2 you would like to move, please open a support ticket with Research Technology and they will assist you.
