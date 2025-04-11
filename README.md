# Tufts University Research Technology Guides

Source repository for the [Tufts University Research Technology Guides][guides-url] -- a collection of resources, documentation, and asynchronous tutorials to advance computational research across disciplines. Developed and maintained by the Research Technology (RT) team within Tufts Technology Services (TTS) at Tufts University.

<!-- prettier-ignore -->

> [!IMPORTANT]
> **See the [published website][guides-url] for content.** The following is a development guide intended for contributors.
> Feedback, suggestions, and error reports can be submitted via <tts-research@tufts.edu> or by [creating an issue][issues-url].

## Table of Contents

- [Prerequisites](#prerequisites)
  - [Miniforge](#miniforge)
  - [GitHub CLI](#github-cli)
  - [Recommended Software](#recommended-software)
- [Initial Setup](#initial-setup)
  - [Optional Pre-Commit Configuration](#optional-pre-commit-configuration)
- [Repository Overview](#repository-overview)
  - [Feature Branches](#feature-branches)
  - [Branch Structure](#branch-structure)
    - [File and Directory Names](#file-and-directory-names)
    - [Placement of Images](#placement-of-images)
    - [Build Artifacts](#build-artifacts)
- [Utility Scripts](#utility-scripts)
  - [`autobuild`](#autobuild)
  - [`build`](#build)
  - [`clean`](#clean)
- [Contribution Workflow](#contribution-workflow)
  - [Creating a Feature Branch](#creating-a-feature-branch)
  - [Committing Changes](#committing-changes)
  - [Submitting a Pull Request](#submitting-a-pull-request)
  - [Merging the Pull Request](#merging-the-pull-request)
- [Publishing Workflow](#publishing-workflow)
  - [Preparing for Publishing](#preparing-for-publishing)
  - [Submitting a Pull Request](#submitting-a-pull-request-1)
  - [Merging the Pull Request](#merging-the-pull-request-1)

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
- `source/_toc.yaml` -- [Sphinx External ToC][toc-url] site map configuration file
- `source/404.md` -- custom 404 page template
- `source/conf.py` -- [Sphinx][sphinx-url] configuration file
- `source/index.md` -- documentation/website root (default landing page)
- `utils` -- various internally-developed utility scripts

All other contents of `source` define the documentation structure and content, with the directory tree corresponding to the site map and files serving as content sources. See corresponding sections below for more information.

#### File and Directory Names

All file and directory names should be URL-friendly and hence only consist of lowercase letters, hyphens (`-`), and numbers. No more than one hyphen may be used in a row and file names should not exceed 32 characters (excluding extension). Prefixes consisting of periods (`.`) or underscores (`_`) are allowed to denote special files and any numerical prefixes should be formatted to two digits. Underscores (`_`) should be used instead of hyphens when naming Python (`py`) scripts. These rules are enforced via the [AutoSlug][autoslug-url] pre-commit hook with automatic fixes applied whenever possible. Any files that do not have an extensions recognized by [AutoSlug][autoslug-url] are treated as directories. File extensions are case-sensitive and care should be take to ensure that the extensions of R scripts (`R`) and R Markdown documents (`Rmd`) are properly capitalized. See below for examples of appropriate names.

```
_special-file.yaml
.hidden-directory
01-file-with-numeric-prefix.rst
python_script.py
python-notebook.ipynb
r-markdown-file.Rmd
this-is-a-directory
this-is-a-file.md
```

#### Placement of Images

Any images intended for inclusion in the source files should be placed in an `img` directory within the same directory as the corresponding source file. Images should be in PNG or SVG format and resized to the desired dimensions whenever possible. All references to images within source files should be relative.

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

- `autobuild` -- automatic self-updating preview build of website hosted on a local server
- `build` -- one-time build of static HTML files
- `clean-autobuild` -- runs `clean` and then `autobuild`
- `clean-build` -- runs `clean` and then `build`
- `clean` -- removal of all build artifacts (needed to ensure a clean build)

New utility scripts should follow the example of existing scripts and be executable without errors from anywhere within the repository. Bash scripts should be developed first and analogous courtesy scripts for Command Prompt and PowerShell users provided when possible. Utility scripts should exit with zero for success and an appropriate positive exit code for failure.

### `autobuild`

The `autobuild` utility script uses [`sphinx-autobuild`][sphinx-autobuild-url] to display an automatically self-updating preview build of the website by running the following command. (The variable `$root` refers to the repository root.)

```sh
sphinx-autobuild --nitpicky --ignore "$root/source/tags" -- "$root/source" "$root/build"
```

The live preview automatically updates whenever chances to source files are detected and can be accessed via [127.0.0.1:8000](http://127.0.0.1:8000) (_localhost_ port number 8000). Note that changes to configuration and template files might not be detected and usually require a clean build to be properly displayed.

### `build`

The `build` utility script uses the standard `sphinx-build` command to build static HTML files that make up the website by running the following command. (The variable `$root` refers to the repository root.)

```sh
sphinx-build --nitpicky "$root/source" "$root/build"
```

The generated HTML files can be previewed by directly opening the landing page (`build/index.html`) or any other desired page using a web browser. Note that some functionality that requires a web server might not work properly when previewing static HTML files.

### `clean`

The `clean` utility script attempts to delete all build artifacts. The usual build process first checks for build artifacts and then only rebuilds the pages where a change to the source file is detected. This means that changes that affect several pages like modifications to the configuration, table of contents, style sheets, or templates might not be accurately reflected when preexisting build artifacts are detected. Hence it is strongly recommended to run `clean` before `autobuild` or `build` whenever making modifications that are not confined to specific source files. Note that the `clean-autobuild` and `clean-build` scripts can be used instead of manually running `clean` before the desired build script.

## Contribution Workflow

The following contribution workflow allows various development efforts to take place simultaneously and ensures changes and new content are incorporated into the `develop` branch in a safe and controlled manner that reduces errors and merge conflicts.

### Creating a Feature Branch

Switch to the `develop` branch and ensure it is up to date.

```
git switch develop
git pull
```

Create a new feature branch from the latest development state.

```
git switch --create <name>
```

Push the new branch to the remote repository and set the upstream tracking branch.

```
git push --set-upstream origin <name>
```

### Committing Changes

Frequent commits are encouraged. The subject line of the commit message should be informative and not exceed500 characters. Commit message bodies are not required but encouraged to provide additional detail where needed. Avoid using sentence-style capitalization and punctuation in the commit message subject line. When using pre-commit hooks, avoid passing the commit message directly via the `-m` flag in case checks fail and automatic fixes are applied.

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

The pull request can be merged when all status checks pass. Run the following to have this automatically happen. Manual intervention is needed when status checks fail and automatic fixes are not possible. Ask for help if needed.

```
gh pr merge --auto
```

The merging of the PR triggers a build and deployment of the development version of the website. Note that this could take several minutes. Once the updated website has been deployed, make sure to navigate to the [development version of the website][guides-dev-url] and sure all changes are reflected as expected.

## Publishing Workflow

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

Thoroughly review all changes and fix any uncovered issues. Request fixes or explanations from original committers when needed.

### Submitting a Pull Request

Create a pull request (PR) to merge the `staging` branch into the `main` branch.

```
gh pr create --base main
```

### Merging the Pull Request

The PR can be merged once it is manually approved by an administrator and all checks pass. Run the following to have this automatically happen. Merging of the PR will trigger the publishing workflow and result in a new build of the public website.

```
gh pr merge --auto
```

## Structure Configuration

Documentation structure is managed using the [Sphinx External ToC][toc-url] extension with the `_toc.yaml` configuration file written such that the site map mimics the layout of the `source` directory. Content is grouped into primary sections with each section appearing in the top navigation bar and having an index file serving as the section root. Primary sections contain content pages which can be further divided into subtrees. Pages in each subtree are ordered using the [natural sort order](https://en.wikipedia.org/wiki/Natural_sort_order) of the source file names. Content pages could also have child pages, in which case their structure resembles that of a primary section with an index file serving as the parent page.

Content pages can be added to preexisting sections, subtrees, and parent pages without having to modify the site map configuration file. Only when adding a new section, subtree, or parent page does the `_toc.yaml` file need to be updated. See the sample `source` directory tree below along with its corresponding site map configuration file for examples on how to define various structures. Note that the `title` field defines how the name of a primary section is displayed in the top navigation bar and the `caption` field defines how the name of a subtree is displayed in the ToC. Content page display names in the secondary sidebar and the ToC are equivalent to their first heading.

```
📂source
 ┣ 📄index
 ┗ 📂primary-section
    ┣ 📄index
    ┣ 📄01-content-page
    ┣ 📄02-content-page
    ┣ 📂10-page-with-children
    ┃  ┣ 📄index
    ┃  ┣ 📄01-child-page
    ┃  ┗ 📄02-child-page
    ┣ 📄21-content-page
    ┣ 📄22-content-page
    ┣ 📄31-subtree-page
    ┗ 📄32-subtree-page
```

```yaml
root: index
subtrees:
  - caption: Primary Section Display Name in ToC
    entries:
      - file: primary-section/index
        title: Primary Section Display Name in Navigation
        subtrees:
          - entries:
              - glob: primary-section/0*
              - file: primary-section/10-page-with-children/index
                entries:
                  - glob: primary-section/10-page-with-children/*
              - glob: primary-section/2*
              - caption: Section Subtree Display Name
                entries:
                  - glob: primary-section/3*
```

File extensions should be omitted when listing source files in the `_toc.yaml` file. This allows for the easy change of source file type without having to modify the structure configuration file. Use file prefixes instead of directories to create subtrees. This avoids the creation of _dead_ URLs where an index file would usually be expected.

[autoslug-url]: https://github.com/TuftsRT/autoslug
[gh-cli-url]: https://cli.github.com/
[guides-dev-url]: https://rtguides.it.tufts.edu/dev/
[guides-url]: https://rtguides.it.tufts.edu/
[issues-url]: https://github.com/TuftsRT/guides/issues
[jinja-url]: https://jinja.palletsprojects.com
[miniforge-url]: https://github.com/conda-forge/miniforge
[myst-md-ext-url]: https://marketplace.visualstudio.com/items?itemName=ExecutableBookProject.myst-highlight
[pre-commit-url]: https://pre-commit.com/
[sphinx-autobuild-url]: https://github.com/sphinx-doc/sphinx-autobuild
[sphinx-url]: https://www.sphinx-doc.org
[toc-url]: https://sphinx-external-toc.readthedocs.io
[vs-code-url]: https://code.visualstudio.com/
