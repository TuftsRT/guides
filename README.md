# Tufts University Research Technology Guides

Collection of resources, documentation, and asynchronous tutorials to advance computational research across disciplines. Developed and maintained by the Research Technology (RT) team within Tufts Technology Services (TTS) at Tufts University.

<!-- prettier-ignore -->

> [!IMPORTANT]
> **See the [published website](https://rtguides.it.tufts.edu) for content.** The following is a development guide intended for staff.

## Prerequisites

### Miniforge

An installation of a Python distribution that supports virtual environments and can handle Conda environment specification YAML files is required. [Miniforge][miniforge-url] is strongly recommended but other similar distributions like Anaconda are acceptable. It is assumed that a conda-enabled console is used to run the example commands in this walkthrough.

### GitHub CLI

The walkthrough also includes commands from the [GitHub CLI][gh-cli-url] (command-line interface), which can be used as a substitute to some `git` commands and browser-based GitHub workflows like creating or merging pull requests (PRs). Once installed, run the following command for initial setup. Authentication via SSH (secure shell) key is strongly recommended over other alternatives.

```
gh auth login
```

### Recommended Software

The majority of the content within the guides is written using MyST (Markedly Structured Text) Markdown, which is an extension of the popular CommonMark Markdown specification and heavily inspired by R Markdown. MyST Markdown contains several special directives that most Markdown editors are unable to properly format or preview. Hence the use of [Visual Studio Code][vs-code-url] along with the [MyST-Markdown][myst-md-ext-url] extension is strongly recommended.

## Initial Setup

Clone the repository using the GitHub CLI or some other equivalent method.

```
gh repo clone tuftsrt/guides
```

Change into the directory containing the cloned repository.

```
cd guides
```

Create the required Python environment using the provided YAML configuration file. (The command provided uses the faster `mamba` package manager included in Miniforge, but the traditional `conda` package manager can be used as well.)

```
mamba env create --file environment.yaml
```

Remember to activate the environment before proceeding. (In certain command line environments, the `mamba` command can also be used for environment activation.)

```
conda activate guides
```

It is assumed that all of the following commands are run from within the `guides` environment. Make sure to activate the environment whenever returning to work with the repository.

### Optional Pre-Commit Configuration

Select [pre-commit][pre-commit-url] hooks are used to ensure consistent style and formatting across all files, fix simple errors like common misspellings, and filter out content that should not be versioned like Jupyter Notebook cell outputs. These hooks are automatically run on all pull requests (PRs) and merging is not possible until all automated checks pass. The PR pre-commit process is able to automatically fix most minor issues that might prevent merging, and hence the use and configuration of pre-commit hooks in the local repository is not required. However, the local use of pre-commit is strongly recommended as the this ensures feature branches are clean and follow a consistent style, and thus are easier to work with. The `pre-commit` library is included in the development environment and can be set up as follows.

Set up all the pre-commit hooks specified in the repository configuration. Run this from the repository root.

```
pre-commit install
```

Test the installation by running the pre-commit hooks against all files in the repository.

```
pre-commit --run --all-files
```

Now the hooks will run before every commit attempt. If any hooks fail or perform automatic fixes, the commit attempt is interrupted and a new commit attempt must be made.

## Repository Overview

The repository contains the following permanent (non-deletable) branches.

- `announcement` -- contains an HTML file used to specify the content of the announcement banner
- `develop` -- receiver of development PRs and the source for the development version of the built website
- `gh-pages` -- receiver of build artifacts and deployment source of the build website
- `main` -- source for the published version of the built website
- `switcher` -- contains a JSON file used to configure the version switcher dropdown
- `staging` -- staging area for content to be pushed to the `main` branch and published

All pushes to the `gh-pages` branch are automated, and developers should **never** interact with this branch directly. Pushes to the `develop` and `main` branches are only allowed via pull requests (PRs), which automatically trigger a new build of the website. The `announcement` and `switcher` branches accept direct pushes and do not trigger a rebuild of the website. However, the content of these branches are directly incorporated into the published website immediately after push. Therefore, any changes to the `announcement` and `switcher` branches should be made with care and the proper functioning of the website manually verified after every push.

The repository may also contain the following non-permanent (deletable) branches.

- `hotfix` -- used for critical changes pushed directly to `main` and then propagated elsewhere
- `pre-commit-ci-update-config` -- automated updates to the pre-commit configuration file

Updates to `pre-commit-ci-update-config` are automated and developers should not directly interact with this branch.

### Feature Branches

All other branches in this repository are _feature branches_ forked from the `develop` branch. These are intended for the active development of a specific feature and all development activity should be confined to these branches. The term _feature_ here could refer to any of the following.

- updates to a specific section of the website
- development of a new section of the website
- changes to configuration files or templates
- development of internal extensions or utility scripts
- updates to the README or other internal documentation

Note that the list above is not exhaustive. Feature branches should be named descriptively (so that it is immediately clear what kind of active development is taking place in the branch) but concisely (not exceeding a soft limit of roughly 30 characters). All branch names should only contain lowercase letters, hyphens (`-`), and numbers.

### Branch Structure

All branches that contain source files (`develop`, `hotfix`, `main`, `staging`, and all feature branches) are structured similarly and contain the following.

- `.github/workflows` -- configuration files for automatic GitHub workflows
- `environment.yaml` -- build environment specification file
- `source/_ext` -- internally-developed and other custom Sphinx extensions
- `source/_static/css` -- cascading style sheet (CSS) files used to override default styling
- `source/_static` -- static HTML content like the website logo and favicon
- `source/_templates` -- custom [Jinja][jinja-url] templates including new templates and default overrides
- `source/_toc.yml` -- [Sphinx External ToC][toc-url] site map configuration file
- `source/404.md` -- custom 404 page template
- `source/conf.py` -- [Sphinx][sphinx-url] configuration file
- `source/index.md` -- documentation/website root (default landing page)
- `utils` -- various internally-developed utility scripts

All other contents of `source` define the documentation structure and content, with the directory tree corresponding to the site map and files serving as content sources. See corresponding sections below for more information.

#### Placement of Images

Any images intended for inclusion in the source files should be placed in an `img` directory within the same directory as the corresponding source file.

#### Build Artifacts

The build process generates the following git-ignored directories that should not be manually modified but are safe to remove.

- `build` -- all build artifacts
- `jupyter_execute` -- executed Jupyter Notebooks derived from source files
- `source/tags` -- automatically generated source files for the tags index

## Utility Scripts

Utility scripts are available for various console environments and can be found in the corresponding directory within the `utils` directory.

- `utils/bash` -- shell scripts (`sh`) prefixed with `#!/bin/bash` and intended to be run via Bash
- `utils/cmd` -- Windows Command Prompt (`cmd`) scripts intended for Miniforge Prompt users
- `utils/pwsh` -- PowerShell (`ps1`) scripts intended for users of cross-platform PowerShell (version 7.X or higher)

All utility scripts are written such that they can be run from anywhere within the repository without errors. (The script uses git commands to derive the repository root and all paths in the script are relative to the repository root.) The following utility scripts are provided for all platforms.

- `autobuild` -- automatic self-updating preview build of website hosted on local server (_localhost_)
- `build` -- one-time build of static HTML files
- `clean-autobuild` -- runs `clean` and then `autobuild`
- `clean-build` -- runs `clean` and then `build`
- `clean` -- removal of all build artifacts (needed to ensure a clean build)

New utility scripts should follow the example of existing scripts and be executable without errors from anywhere within the repository. Bash scripts should be developed first and analogous courtesy scripts for Command Prompt and PowerShell users provided when possible.

## Contribution

The following contribution workflow allows various development efforts to take place simultaneously and ensures changes and new content are incorporated into the `develop` branch in a safe and controlled manner that reduces errors and merge conflicts.

### Creating a Feature Branch

Switch to the `develop` branch and ensure it is up to date.

```
git switch develop
git pull
```

Create a new feature branch from the latest development state.

```
git switch --create [your feature branch name]
```

Push the new branch to the remote repository and set the upstream tracking branch.

```
git push --set-upstream origin [your feature branch name]
```

### Submitting a Pull Request

Create a pull request (PR) to merge your feature branch into the `develop` branch. Follow the prompts to add a title and description for your PR.

```
gh pr create --base develop
```

View the status and details of your PR.

```
gh pr view
```

You can also open the PR on the GitHub website.

```
gh pr view --web
```

### Merging the Pull Request

The pull request can be merged when all status checks pass. Run the following to have this automatically happen. Manual intervention is needed when status checks fail. Ask for help if needed.

```
gh pr merge --auto
```

## Publishing

### Preparing for Publishing

Switch to the `staging` branch to prepare for publishing.

```
git switch staging
git pull
```

Merge the latest version `develop` branch into `staging` to incorporate the changes.

```
git fetch origin
git merge origin/develop
```

### Creating a Pull Request

Create a pull request (PR) to merge the `staging` branch into the `main` branch.

```
gh pr create --base main
```

### Merging the Pull Request

The PR can be merged once all checks pass and it is manually approved by an administrator and all checks pass. Run the following to have this automatically happen. Merging of the PR will trigger the publishing workflow and result in a new build of the public website.

```
gh pr merge --auto
```

[gh-cli-url]: https://cli.github.com/
[jinja-url]: https://jinja.palletsprojects.com
[miniforge-url]: https://github.com/conda-forge/miniforge
[myst-md-ext-url]: https://marketplace.visualstudio.com/items?itemName=ExecutableBookProject.myst-highlight
[pre-commit-url]: https://pre-commit.com/
[sphinx-url]: https://www.sphinx-doc.org
[toc-url]: https://sphinx-external-toc.readthedocs.io
[vs-code-url]: https://code.visualstudio.com/
