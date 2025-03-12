# Globus

> This step-by-step guide will show you how to log into Globus and use it to transfer files reliably and securely. You will become familiar with basic Globus concepts and interfaces, and begin to experience how Globus can help you spend more time on your research and less time on data management.

## What is Globus?

[Globus](https://www.globus.org/) is research cyberinfrastructure, developed and operated as a not-for-profit service by the University of Chicago.

Globus provides a secure, unified interface to your research data. Use Globus to 'fire and forget' high-performance data transfers between systems within and across organizations.

## Key Concept: _endpoint_

An endpoint is a server that hosts collections. If you want to be able to access, share, transfer, or manage data using Globus, the first step is to create an endpoint on the system where the data is (or will be) stored.

[Globus Connect](https://www.globus.org/globus-connect) is used to create endpoints. An endpoint can be a laptop, a personal desktop system, a laboratory server, a campus data storage service, a cloud service, or an HPC cluster. As explained below, it’s easy to set up your own Globus endpoint on a laptop or other personal system using Globus Connect Personal. Administrators of shared services (like campus storage servers) can set up multi-user endpoints using Globus Connect Server. You can use endpoints set up by others as long as you’re authorized by the endpoint administrator or by a collection manager.

### Globus File Transfer Prerequisites

> Globus Account(s)
>
> Two Endpoints
>
> VPN or Tufts network is not required

```{gallery-grid}
---
grid-columns: 1
---
- header: "{fas}`book` Globus Account Setup"
  content: "Tufts has a subscription to Globus, and you can set up a Globus account with your Tufts credentials. You can also link other accounts, either personal or through other institutions. Link Tufts Account"
  link: "42-globus-account-setup.html"

- header: "{fas}`book` Tufts Globus Collections"
  content: "Collections create by Tufts Research Technology, allows you to share and transfer files to and from your storage on Tufts HPC Cluster and Sharepoint."
  link: "43-globus-tufts-collections.html"

- header: "{fas}`book` Local Endpoint - Globus Connect Personal"
  content: "Create a private collection on your own computer. Globus Connect Personal allows you to share and transfer files to and from your Mac/Linux/Windows laptop or desktop computer."
  link: "44-globus-connect-personal.html"

- header: "{fas}`book` Transfer Files between Two Endpoints"
  content: "How to start and monitory your file transfers."
  link: "45-globus-transfer-files.html"

```

[Learn More](https://docs.globus.org/guides/tutorials/manage-files/transfer-files/)
