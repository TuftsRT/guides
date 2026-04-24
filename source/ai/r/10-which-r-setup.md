---
tags: r data-science rstudio
---

# Which R Setup is Right for You?

There are several ways to access and use R, and choosing the right one can make your work significantly easier. This guide outlines the options we recommend for Tufts researchers and students, organized by experience level and research needs.

## Our Recommendations

- **Beginners and first-time R users:**
  - [Posit Cloud](#posit-cloud-a-beginner-friendly-browser-based-option): browser-based, no installation required
- **Intermediate to advanced users and most research workflows:**
  - [RStudio Desktop with renv](#rstudio-desktop-our-recommended-local-setup): local installation, full control over your environment
- **Users with large datasets or computationally intensive analyses:**
  - [R on the Tufts HPC Cluster](#r-on-the-tufts-hpc-cluster): access to high-performance computing resources

Read the sections below to understand the pros and cons of each option, then follow the links to the relevant setup instructions.

---

## Posit Cloud: A Beginner-Friendly, Browser-Based Option

[Posit Cloud](https://posit.cloud/) (formerly RStudio Cloud) is a browser-based version of the RStudio interface. It lets you use R and RStudio immediately, without installing anything on your computer.

Posit Cloud is a good choice if:

- You are new to R and want to start quickly without worrying about installation
- You are working on a managed or restricted computer where software installation is difficult
- You want to share a project with collaborators and have everyone work in the same environment

**Limitations:**

- Free accounts have monthly compute-hour limits, which can be reached quickly for intensive analyses
- Not suitable for sensitive or IRB-restricted data (data is uploaded to cloud servers)
- Less control over your R and package versions compared to a local installation

> **Data privacy:** Do not use Posit Cloud with sensitive, IRB-restricted, or HIPAA-regulated data. Consult the [Tufts data storage finder](https://access.tufts.edu/data-finder) or contact datalab-support@tufts.edu if you are unsure about where to store your data.

### Summary of Pros and Cons for Posit Cloud

**Pros:**

- No installation required; works in any web browser
- Familiar RStudio interface
- Easy project sharing and collaboration

**Cons:**

- Monthly compute limits on free tier
- Not appropriate for sensitive or IRB-restricted data
- Less control over package versions and environment management
- Requires an internet connection

---

## RStudio Desktop: Our Recommended Local Setup

For most research workflows, we recommend installing **R** and **RStudio Desktop** on your local computer, combined with the **renv** package for managing your project libraries.

**R** is the underlying programming language and statistical engine. **RStudio Desktop** is an IDE (Integrated Development Environment) that provides a user-friendly interface for writing R code, viewing outputs and plots, and managing your files and environment. **renv** is a package for creating isolated, reproducible package environments within individual R projects, similar in purpose to conda environments in Python.

This setup is the standard for R users in academia and industry. It gives you full control over your R and package versions, works well with sensitive data (everything stays local), and makes your analyses reproducible.

See [RStudio Desktop Setup](20-rstudio-desktop.md) for installation and getting-started instructions.

### Summary of Pros and Cons for RStudio Desktop

**Pros:**

- Free for all users and use cases
- Cross-platform (Windows, macOS, Linux)
- Works with sensitive or IRB-restricted data (data stays on your computer)
- Full control over R and package versions
- Reproducible project environments with renv
- Large user community and extensive documentation
- Excellent support for R Markdown and Quarto for reproducible reports

**Cons:**

- Requires installation and some initial setup
- Not suitable for very large datasets or computationally intensive jobs (limited by your local hardware)

---

## R on the Tufts HPC Cluster

If your analyses require more memory or computing power than your personal computer can provide, such as large genomic datasets, extensive simulation studies, or long-running statistical models, you may want to use R on the Tufts High-Performance Computing (HPC) cluster.

The HPC gives you access to large amounts of memory and CPU cores. There are two main ways to use R on the HPC:

- **RStudio Server via OnDemand**: a browser-based RStudio interface running on the cluster, good for interactive analyses
- **Batch jobs via SLURM**: submit R scripts to run non-interactively, good for long-running or parallelized analyses

For full instructions, see [R on the Tufts HPC](30-r-on-hpc.md).

### Summary of Pros and Cons for R on the HPC

**Pros:**

- Access to large memory and many CPU cores
- Suitable for computationally intensive analyses
- RStudio Server provides a familiar interface
- SLURM batch jobs allow long-running or automated analyses

**Cons:**

- Requires an HPC account (apply at [https://it.tufts.edu/high-performance-computing](https://it.tufts.edu/high-performance-computing))
- Some additional complexity in job submission and file management
- Interactive sessions have resource limits and time caps
- Not appropriate for all data types; check HPC data policies for sensitive data
