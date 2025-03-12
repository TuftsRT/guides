---
tags: bioinformatics
---

# Using containers with Singularity to ensure reproducible research

Shirley Li, Bioinformatician, TTS Research Technology
xue.li37@tufts.edu

Date: 2024-11-01

**We recommend you read this paper before delving into the details: [Containerized Bioinformatics Ecosystem for HPC](https://ieeexplore.ieee.org/abstract/document/10027523)**

Singularity is pivotal in bioinformatics for ensuring computational reproducibility and facilitating the seamless transition of workflows from development to production environments. This document delves deeper into specific applications of Singularity, focusing on BioContainers and its utility within Nextflow pipelines.

## BioContainers with Singularity

[BioContainers](https://biocontainers.pro/) is a project that provides pre-built container images for bioinformatics software, facilitating easy deployment of thousands of software packages. Singularity is often used to run these containers because it combines compatibility with security in multi-user environments like shared computing resources.

### Key Benefits:

- Standardization: BioContainers standardizes the deployment of bioinformatics software, which significantly reduces inconsistencies that arise from software installation issues.
- Accessibility: Researchers can access and run any software provided by BioContainers without concerns about dependencies or conflicts with other installed software.

## Integration of Singularity in Nextflow Pipelines

[Nextflow](https://www.nextflow.io/) is a workflow management system that simplifies the parallelization and execution of bioinformatics pipelines across multiple compute infrastructures. Singularity is a preferred containerization technology used in Nextflow due to its robustness and compatibility with HPC environments.

### Workflow:

- Pipeline Configuration: Users specify Singularity as the container engine in the Nextflow configuration file.
- Container Retrieval: Nextflow automatically pulls the required Singularity container from BioContainers or other container registries.
- Execution: The pipeline tasks are executed within the containers, ensuring that each step of the workflow runs within a controlled and consistent environment.

## Reference

- Singularity Documentation: Visit the [Singularity](https://sylabs.io/singularity/) website for comprehensive guides and documentation.
- BioContainers Repository: Access thousands of bioinformatics containers at [BioContainers](https://biocontainers.pro/).
- Nextflow and Singularity: Learn how to integrate Singularity in Nextflow workflows from the [Nextflow documentation](https://www.nextflow.io/docs/latest/singularity.html).
- [nf-core pipeline configuration](https://nf-co.re/docs/usage/getting_started/configuration) with Singularity and other container engines.

## Previous workshop

[Pulling, running and building containerized HPC applications](https://zhan4429.github.io/TuftsContainers.github.io/)
