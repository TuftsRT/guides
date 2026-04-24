---
tags: r rstudio renv
---

# RStudio Desktop: Installation and Setup

This guide walks you through installing R and RStudio Desktop on your local computer and getting your first project set up. It also introduces **renv**, our recommended tool for managing R package environments.

## Step 1: Install R

R and RStudio are two separate pieces of software. You need to install R first, then RStudio.

1. Go to [https://cran.r-project.org/](https://cran.r-project.org/)
1. Click the link for your operating system (Windows, macOS, or Linux)
1. Download and run the latest installer
   - **Windows:** Click "base", then download and run the `.exe` installer
   - **macOS:** Download the `.pkg` file appropriate for your chip (Apple Silicon M1/M2/M3 or Intel)
   - **Linux:** Follow the instructions for your distribution

> **Note:** The Tufts Data Lab also maintains instructions for installing R and RStudio at [https://sites.tufts.edu/datalab/installing-r-and-rstudio/](https://sites.tufts.edu/datalab/installing-r-and-rstudio/).

## Step 2: Install RStudio Desktop

RStudio Desktop is a free IDE (Integrated Development Environment) that provides a user-friendly interface for writing and running R code.

1. Go to [https://posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/)
1. Download the free **RStudio Desktop** installer for your operating system
1. Run the installer

Once installed, open RStudio. It will automatically detect your R installation.

## The RStudio Interface

When you open RStudio, you will see four panes:

- **Source** (top left): Write and edit your R scripts and R Markdown files
- **Console** (bottom left): Run R commands interactively and see output
- **Environment/History** (top right): View loaded variables and command history
- **Files/Plots/Packages/Help** (bottom right): Browse files, view plots, manage packages, and access documentation

## Step 3: Create an RStudio Project

RStudio Projects are a core organizing concept: each project gets its own working directory, history, and (with renv) its own package environment. We strongly recommend using a project for each analysis or research task.

1. In RStudio, go to **File → New Project**
1. Choose **New Directory → New Project**
1. Give your project a name and choose where to save it
1. Click **Create Project**

Your project will open in a fresh RStudio session. RStudio will automatically set the working directory to your project folder, so file paths in your code can be relative to that folder.

## Step 4: Install and Use Packages

R's power comes largely from its ecosystem of packages. To install a package, use `install.packages()` in the Console:

```r
install.packages("tidyverse")
```

To load a package and use it in your session:

```r
library(tidyverse)
```

You only need to install a package once (per R version), but you need to call `library()` at the start of each session or script that uses it.

## Step 5: Manage Package Environments with renv

As you work on research projects over time, you may find that updates to packages cause your old code to break. **renv** solves this by creating a self-contained, project-specific package library. When you share your project or return to it later, renv can restore the exact package versions you used.

### Initialize renv in Your Project

In your project's Console, run:

```r
install.packages("renv")
renv::init()
```

This creates a project-local library and a `renv.lock` file that records all your package versions. From this point on, packages you install while in the project will be stored in the project library, not your global R library.

### Save Your Environment

After installing or updating packages, save the current state of your environment:

```r
renv::snapshot()
```

This updates `renv.lock` to reflect the current package versions.

### Restore an Environment

To restore a saved environment (for example, when reopening the project on a new computer or sharing with a collaborator):

```r
renv::restore()
```

renv will install the exact package versions recorded in `renv.lock`.

> **Tip:** Commit `renv.lock` to version control (e.g., Git) along with your code. This allows collaborators and your future self to reproduce your exact environment.

## Writing and Running Code

### R Scripts (.R files)

R scripts are plain text files containing R code. To create one:

1. Go to **File → New File → R Script**
1. Write your code in the Source pane
1. Run a single line with **Ctrl+Enter** (Windows/Linux) or **Cmd+Enter** (macOS)
1. Run the entire script with **Ctrl+Shift+Enter** / **Cmd+Shift+Enter**

### R Markdown (.Rmd files) and Quarto (.qmd files)

R Markdown and Quarto let you combine code, text, and output in a single document, making it well suited for reproducible reports and analyses. To create an R Markdown document:

1. Go to **File → New File → R Markdown**
1. Choose an output format (HTML, PDF, Word)
1. Click **Knit** to run the code and render the document

Quarto (`.qmd`) is the next-generation version of R Markdown and supports both R and Python. RStudio supports Quarto natively.

## Next Steps

Once you have R, RStudio, and renv set up, you are ready to start your analysis. Browse the tutorial pages in this section for hands-on guides covering data visualization, statistical modeling, and more.

For help, email the Tufts Data Lab at **datalab-support@elist.tufts.edu**.
