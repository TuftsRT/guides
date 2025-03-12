---
tags: bioinformatics
---

# Running nfcore pipelines on Tufts HPC

Shirley Li, Bioinformatician, TTS Research Technology
xue.li37@tufts.edu

Date: 2024-11-01

You have two options for running the nf-core pipeline on Tufts HPC: through Open OnDemand or the command-line interface.

## 1. Open Ondemand

Navigate to [Open Ondemand](https://ondemand.pax.tufts.edu/), under `Bioinformatics Apps`, choose the appropriate application for the pipeline you wish to execute.

## 2. Linux Command Line Interface

By using `module avail`, you can see the available modules. There are two ways to run the nf-core pipeline via the command-line interface:

### Option 1

Load the necessary modules

```
module load singularity/3.8.4
module load nf-core/2.13.1
```

Run the pipeline with this command:

```
nextflow run nf-core/mag -profile tufts ...
```

For a comprehensive list of available pipelines, visit the [nf-core website](https://nf-co.re/pipelines). There are currently 107 pipelines available through nf-core.

### Option 2

Load singularity module and the specific pipeline module:

```
module load singularity/3.8.4
module load nf-core-mag/2.5.4
```

Execute the pipeline using:

```
mag ...
```

Currently, 31 pipelines have been deployed as modules on Tufts HPC for ease of use.
