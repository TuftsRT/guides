# HPC Policy

## Users

The HPC cluster is available for use by all Tufts Faculty, Staff and Students.
Faculty with a Tufts University appointment can request additional storage for use by their lab.

## Acceptable Use Policy

Tufts HPC cluster is an "institutional systems" and its acceptable use is governed by the university wide ["Use of Information Systems Policy"](https://it.tufts.edu/sites/default/files/uploaded-files/2018-09/use-institutional_systems_2018-06-18_0.pdf) policy.

## Data

- Data storage on the HPC Cluster is provided by the Tier 1 HPC research storage system. Details on this service can be found at (https://it.tufts.edu/research-data-storage).
- No restricted data is allowed on Tufts HPC cluster.

### Scratch Space
Scratch disk space is a commonly implemented feature in high performance computing environments. This temporary space allows users the freedom needed to work with large spikes in data usage typical of HPC. It promotes storing temporary and transient data such as check points, uncompressed copies of data sets, and caches in a space specifically designed for it. Data on scratch is automatically deleted if not accessed in the last 21 days.  

**Benefits of Scratch Space**
* Promotes good data management practices by providing a place that "temporary" data should go that is separate and can be managed differently. 
* Allows uniform documentation and configuration of applications and workflows to utilize this temporary space for transient files. 
* Promotes users keeping transient, and non-important data separate from other data. This allows the scratch space to not have snapshots, backups, or indexing, providing significant costs saving. 
* Allows users without a designated project folder, such as undergraduate students, to complete homework or projects on the HPC system that would otherwise require too much disk space.

**HPC Users**
* Each HPC user has a personal scratch folder available to use
* It is mounted at /cluster/scratch/utln
* The quota is 100GB
* Files in scratch are automatically deleted if not accessed in the last 21 days 

**Researchers**
* Tufts-employed faculty members can request the scratch quota of a given lab member be raised up to 15TB. 
* The allocation of this extra space is subject to review and approval by Research Technology. 

**Implementation Details** 
* Data is automatically deleted if not accessed in 21 days 
* There are no backups or snapshots of data stored in scratch 
* A total of 300 TB of Tier 1 storage has been designated for use as scratch space. In the event this fills up we may revisit the amount of space provided via this policy. 
* In the unlikely event any single user needs a scratch space over 15TB the additional storage allocated will be charged at the Tier 1 storage rate.
* Running `touch` commands or similar operations to modify timestamps and bypass this cleanup policy is prohibited. Users who engage in this behavior will lose the privilege of using scratch storage.

## HPC Researcher Contribution Node

Faculty provided contribute nodes are configured to give the contributing researcher priority access but are otherwise
operated the same as other nodes in the cluster.

- The owners of the equipment are given priority access while any excess capacity is returned to the pool (“preempt partition”) for general, shared use by the Tufts community.
- The owners of the equipment determine the access and queuing policies for their portion of the HPC.
- When not in use, contrib nodes provide compute resources to the general Tufts community via a preempt partition.
- When node owners require access to their contrib nodes, non-owner jobs are preempted within approx. 30 seconds to ensure owners have timely access to contributed capacity.

Node types currently available for purchase can be found in our [Contribute Node Catalog](contribute-nodes.md).

### Lifecycle

As with all technology, HPC nodes have a service life. As the compute nodes age, they lose vendor and operating system support, creating maintenance difficulties and security risks. To ensure that Tufts researchers receive sufficient advance notice for the planning of HPC compute node end-of-life, TTS will keep researchers informed of the state of the lifecycle of their contributed nodes, allowing a total service life of the warranty period and a best-effort support period lasting a maximum of 7 years. At purchase, every HPC compute node must be procured with an associated warranty. Nodes that fail without warranty support will be decommissioned; and, to ensure the best value from contributed nodes, TTS recommends a 5-year warranty. This policy is intended to help researchers avoid significant unplanned interruptions to their research caused by system failures. TTS Research Technology will keep researchers informed about the state of their nodes and work with researchers to determine node retirement dates and replacement estimates.

### Installation Costs

Effective January 1, 2025, installation costs for researcher-purchased compute nodes at the Massachusetts Green High Performance Computing Center (MGHPCC) will be funded directly by the school/department or faculty purchasing the equipment. Grant proposals submitted and startup packages committed before September 1st, 2024 that include node purchases will be exempted from this change. While costs may vary based on the number of nodes installed, TTS recommends a planning estimate of \$1,000 total per typical installation of 10-20 nodes.

## Acknowledgement

To acknowledge your use of the clusters, we request that you use our [DOI](http://dx.doi.org/10.25833/235d-b074) or the following wording:

The authors acknowledge the Tufts University High Performance Compute Cluster (https://it.tufts.edu/high-performance-computing) which was utilized for the research reported in this paper.

We would also appreciate your informing us of publications in which you acknowledge one of the clusters.
