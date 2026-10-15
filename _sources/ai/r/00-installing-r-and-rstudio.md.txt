---
tags: r rstudio installation
---

# Installing R and RStudio for the First Time

This guide covers setting up R on a personal computer for the first time. The process takes about 15 minutes.

## Why Two Downloads Are Needed

Two separate programs are required:

| Program             | Description                                                                                        | Download Source                                   |
| ------------------- | -------------------------------------------------------------------------------------------------- | ------------------------------------------------- |
| **R**               | The programming language itself, which performs the actual computing.                              | [The R Project][r-project] (via CRAN)             |
| **RStudio Desktop** | An integrated development environment (IDE) for writing and running R code. It requires R to work. | [Posit's RStudio download page][rstudio-download] |

R is the engine and RStudio is the dashboard. RStudio cannot run any code until R is installed. Both programs are free and open source.

```{note}
Posit (the company that develops RStudio) has redesigned its website, and the RStudio page no longer links to R. R must now be downloaded separately from the R Project website. Install **R first**, then RStudio.
```

## Step 1: Install R

1. Go to [r-project.org][r-project].
1. Click **download R** (or **CRAN** in the left-hand menu). CRAN (Comprehensive R Archive Network) is the official network of servers that distributes R.
1. A list of _mirrors_ (copies of the same files hosted around the world) appears. Select **0-Cloud** at the top of the list, which automatically uses a nearby server.
1. Under **Download and Install R**, select the link for the appropriate operating system and follow the corresponding instructions below.

```{tip}
Steps 1–3 can be skipped by going directly to [cloud.r-project.org][cran-cloud].
```

### Windows

1. Click **Download R for Windows**.
1. Click **base** (or **install R for the first time**).
1. Click **Download R-4.x.x for Windows** (currently `R-4.6.1-win.exe`).
1. Open the downloaded `.exe` file and click through the installer. The default options are appropriate for most users.

R requires Windows 10 or later.

### macOS

1. Click **Download R for macOS**.
1. Choose the installer that matches the Mac's processor. To check the processor type, open the Apple menu and select **About This Mac**.
   - **Chip: Apple M1, M2, M3, M4, etc.** (Apple silicon): download the **`arm64`** file (currently `R-4.6.1-arm64.pkg`). Requires macOS 14 (Sonoma) or newer.
   - **Processor: Intel**: download the **`x86_64`** file (currently `R-4.6.1-x86_64.pkg`). Requires macOS 11 (Big Sur) or newer.
1. Open the downloaded `.pkg` file and follow the installer. The default options are appropriate for most users.
   - If macOS requests permission for the installer to access the Downloads folder, click **Allow**.

```{note}
Some R graphics and packages require [XQuartz][xquartz]. Install it if an error about X11 appears later.
```

### Linux

Click **Download R for Linux** and select the appropriate distribution (Ubuntu, Debian, Fedora/Red Hat, etc.). Each distribution has its own instructions, which usually involve installing R through the system package manager (for example, `apt` on Ubuntu). Follow the steps on the CRAN page for the relevant distribution.

## Step 2: Install RStudio Desktop

1. Go to the [RStudio product page][rstudio-product] and click **Download RStudio**, or go directly to the [RStudio download list][rstudio-download].

1. Find the installer for the appropriate operating system:

   | Operating System                         | File to Download                   |
   | ---------------------------------------- | ---------------------------------- |
   | Windows 11                               | `RStudio-2026.xx.x-xxx.exe`        |
   | macOS 14+                                | `RStudio-2026.xx.x-xxx.dmg`        |
   | Ubuntu / Debian                          | `rstudio-2026.xx.x-xxx-amd64.deb`  |
   | RHEL (Red Hat Enterprise Linux) / Fedora | `rstudio-2026.xx.x-xxx-x86_64.rpm` |

1. Install RStudio:

   - **Windows:** Run the `.exe` file and accept the defaults.
   - **macOS:** Open the `.dmg` file and **drag the RStudio icon into the Applications folder**. Do not run RStudio directly from the disk image.
   - **Linux:** Open the `.deb` or `.rpm` file with the system software installer, or install it from the terminal.

```{note}
Only the free **RStudio Desktop** (open source edition) is needed. The **RStudio Desktop Pro** and **Try Posit AI** options are not required for coursework or research.
```

## Step 3: Verify the Installation

1. Open **RStudio**. Opening R separately is not necessary.

1. RStudio displays a window with several panes. In the **Console** pane (usually bottom-left), the first lines should read something like the following.

   ```text
   R version 4.6.1 (2026-06-24) -- "Happy Hop"
   ```

1. Click in the Console next to the `>` prompt, type the following, and press **Enter**.

   ```r
   1 + 1
   ```

   The following output should appear.

   ```text
   [1] 2
   ```

1. Install a package. Packages are add-ons that extend the functionality of R. Type the following.

   ```r
   install.packages("tidyverse")
   ```

   Installation may take a few minutes. When it finishes, load the package.

   ```r
   library(tidyverse)
   ```

   If no error messages appear, R and RStudio are installed and working correctly.

## Troubleshooting

**RStudio cannot find R or prompts for an R installation.**
Either R is not installed or it was installed while RStudio was open. Close RStudio, confirm that Step 1 was completed, and reopen RStudio. On Windows, the R version can be selected under **Tools → Global Options → General → R version**.

**The installer requests an administrator password.**
Tufts-managed computers may not grant administrator rights. Contact Tufts Technology Services (TTS) for help installing the software.

**`install.packages()` asks whether to use a personal library instead.**
Answer **yes**. R saves packages to a folder in the user account, which is normal.

**On Windows, a package fails to install with a message about Rtools or compiling from source.**
Most packages do not require Rtools. If one does, install [Rtools], choosing the version that matches the installed R version (for example, Rtools 4.6 for R 4.6.x).

**On macOS, an error mentions X11 or `tcltk`.**
Install [XQuartz], then restart the computer. Reinstall XQuartz after every major macOS upgrade.

**The wrong macOS version of R (Intel or Apple silicon) was installed.**
Download and run the correct installer. It replaces the previous installation.

## Keeping R and RStudio Up to Date

- **R** releases new versions several times a year. To update, download and install the newest version from [CRAN][cran-cloud] following the same process as Step 1. Reinstall packages afterward, because packages are tied to the major and minor version of R (for example, 4.6).
- **RStudio** displays a notification when an update is available. Updates can also be checked manually under **Help → Check for Updates**.

## Additional Resources

- [R for Data Science (2e)](https://r4ds.hadley.nz/) -- a free online book and a recommended starting point
- [RStudio IDE User Guide](https://docs.posit.co/ide/user/)
- [CRAN R FAQ](https://cran.r-project.org/faqs.html)
- [RT Workshops](https://go.tufts.edu/workshops)
- [cran-cloud](https://cloud.r-project.org/)
- [r-project](https://www.r-project.org/)
- [rstudio-download](https://docs.posit.co/ide/user/#rstudio-ide-oss-downloads)
- [rstudio-product](https://posit.co/products/open-source/rstudio)
- [rtools](https://cloud.r-project.org/bin/windows/Rtools/)
- [xquartz](https://www.xquartz.org/)
