# TTS Research Technology Guides

Source repository for the [TTS Research Technology Guides][guides-url] -- a collection of resources, documentation, and asynchronous tutorials to advance computational research across disciplines. Developed and maintained by the Research Technology (RT) team within Tufts Technology Services (TTS) at Tufts University.

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
- [Structure Configuration](#structure-configuration)
- [Content Types](#content-types)
  - [Narrative](#narrative)
    - [MyST Markdown](#myst-markdown)
    - [reStructuredText](#restructuredtext)
  - [Executable](#executable)
    - [Jupyter Notebook](#jupyter-notebook)
    - [MyST Markdown Notebook](#myst-markdown-notebook)
    - [R Markdown Notebook](#r-markdown-notebook)
- [Style Guidelines](#style-guidelines)
  - [Prose](#prose)
  - [Code](#code)
- [Accessibility](#accessibility)
- [Subject Tags](#subject-tags)
  - [MyST and R Markdown](#myst-and-r-markdown)
  - [reStructuredText](#restructuredtext-1)
  - [Jupyter Notebook](#jupyter-notebook-1)
- [Contribution Workflow](#contribution-workflow)
  - [Creating a Feature Branch](#creating-a-feature-branch)
  - [Updating the Environment](#updating-the-environment)
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

The majority of the content within the guides is written using MyST (Markedly Structured Text) Markdown, which is an extension of the popular CommonMark Markdown specification and heavily inspired by R Markdown. MyST Markdown contains several special directives that most Markdown editors are unable to properly format or preview. Hence the use of [Visual Studio Code][vs-code-url] along with the [MyST-Markdown][myst-md-ext-url] extension is strongly recommended. Alternatively, the JupyterLab distribution included in the project environment is equipped with an extension that adds MyST Markdown support and can be used to edit and preview MyST markdown files instead.

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
mamba env create --file environment.yaml --yes
```

Remember to activate the environment before proceeding. (In certain command line interfaces, the `mamba` command can also be used for environment activation.)

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
pre-commit run --all-files
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

Note that the list above is not exhaustive. Feature branches should be named descriptively (so that it is immediately clear what kind of active development is taking place in the branch) but concisely (not exceeding a soft limit of roughly 30 characters). All branch names should only contain lowercase letters, hyphens (`-`), and numbers. Do not use more than one hyphen in a row.

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

All other contents of `source` define the documentation structure and content, with the directory tree corresponding to the site map and files serving as content sources. See corresponding sections below for more information. Note that all files in `source` that are not listed above are automatically published online even if not linked to from anywhere in the content. Do not place any non-content files like utility scripts or developer-facing documentation within the `source` directory.

#### File and Directory Names

All file and directory names should be URL-friendly and hence only consist of lowercase letters, hyphens (`-`), and numbers. No more than one hyphen may be used in a row and file names should not exceed 32 characters (excluding extension). Prefixes consisting of periods (`.`) or underscores (`_`) are allowed to denote special files and any numerical prefixes should be formatted to two digits. Underscores (`_`) should be used instead of hyphens when naming Python (`py`) scripts. These rules are enforced via the [AutoSlug][autoslug-url] pre-commit hook with automatic fixes applied whenever possible. Any files that do not have an extension recognized by [AutoSlug][autoslug-url] are treated as directories. File extensions are case-sensitive and care should be taken to ensure that the extensions of R scripts (`R`) and R Markdown documents (`Rmd`) are properly capitalized. See below for examples of appropriate names.

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

Any images intended for inclusion in the source files should be placed in an `img` directory within the same directory as the corresponding source file. Images should be in PNG or SVG format and resized to the desired dimensions whenever possible. All references to images within source files should be relative. Do not include references to any externally hosted images

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

## Content Types

Content that is not code-heavy is considered _narrative_ and can be written using either MyST Markdown or reStructuredText. Narrative content can contain code snippets but these are not executed and their output is not intended to be included. Content that is focused on code and intends to include both the code itself and its output is considered _executable_ and can be written using a variety of notebook formats. Code snippets included in executable content are executed during the build and the outputs are automatically included in the built document. Efforts are underway to make the code included in the pages derived from executable content interactive using WebAssembly (WASM).

All content regardless of format must adhere to the following rules.

- **Files must have a single title.** Generally this means that files must begin with a single first-level heading and no other first-level headings should appear. (R Markdown is an exception. There the title is defined in the YAML header and hence multiple first-level headings can be used.)
- **Headings must increase linearly.** Subsections of sections with a first-level header should have a second-level header, the subsections of which should have a third-level header and so on. Skipping a header level will result in a build error.

### Narrative

Narrative content can be written using either MyST Markdown (stored in `md` files) or reStructuredText (stored in `rst` files). The use of MyST Markdown over reStructuredText is strongly encouraged due to its superior functionality and readability. Support for reStructuredText exists primarily to allow the inclusion of preexisting legacy materials.

#### MyST Markdown

MyST (Markedly Structured Text) Markdown is an extension of the popular [CommonMark](https://commonmark.org/) Markdown specification and is heavily inspired by [R Markdown](https://rmarkdown.rstudio.com/). It supports simple formatting, lists, images, tables, mathematical formulas, and code snippets like most other Markdown flavours and adds various roles and directives that allow for extra functionality like admonitions, footnotes, citations, and glossaries. Furthermore, the [PyData Sphinx Theme](https://pydata-sphinx-theme.readthedocs.io/en/stable/index.html) and the [Sphinx Design Extension](https://sphinx-design.readthedocs.io/en/latest/) add additional design elements and functionality like grids, cards, dropdowns, tabs, sidebars, and iconography. See below for relevant resources.

- [CommonMark Markdown Reference](https://commonmark.org/help/)
- [CommonMark Markdown Tutorial](https://commonmark.org/help/tutorial/)
- [MyST Markdown Cheat Sheet](https://jupyterbook.org/en/stable/reference/cheatsheet.html)
- [MyST Markdown Documentation](https://myst-parser.readthedocs.io/en/latest/index.html)
  - [Typography](https://myst-parser.readthedocs.io/en/latest/syntax/typography.html)
  - [Admonitions](https://myst-parser.readthedocs.io/en/latest/syntax/admonitions.html)
  - [Images](https://myst-parser.readthedocs.io/en/latest/syntax/images_and_figures.html)
  - [Tables](https://myst-parser.readthedocs.io/en/latest/syntax/tables.html)
  - [Code](https://myst-parser.readthedocs.io/en/latest/syntax/code_and_apis.html)
  - [Cross-References](https://myst-parser.readthedocs.io/en/latest/syntax/cross-referencing.html)
  - [Math and Equations](https://myst-parser.readthedocs.io/en/latest/syntax/math.html)
- [Sphinx Design Extension Documentation](https://sphinx-design.readthedocs.io/en/latest/)
  - [Grids](https://sphinx-design.readthedocs.io/en/latest/grids.html)
  - [Cards](https://sphinx-design.readthedocs.io/en/latest/cards.html)
  - [Dropdowns](https://sphinx-design.readthedocs.io/en/latest/dropdowns.html)
  - [Tabs](https://sphinx-design.readthedocs.io/en/latest/tabs.html)
  - [Iconography](https://sphinx-design.readthedocs.io/en/latest/badges_buttons.html)
- [Elements Specific to the PyData Sphinx Theme](https://pydata-sphinx-theme.readthedocs.io/en/stable/user_guide/theme-elements.htm)
- [PyData Sphinx Theme Kitchen Sink (preview of almost all design elements)](https://pydata-sphinx-theme.readthedocs.io/en/stable/examples/kitchen-sink/index.html)

MyST Markdown is extremely configurable and has various syntax extensions. The following have been enabled for this project.

- [Typography](https://myst-parser.readthedocs.io/en/latest/syntax/optional.html#typography)
- [Strikethrough](https://myst-parser.readthedocs.io/en/latest/syntax/optional.html#strikethrough)
- [Math Shortcuts](https://myst-parser.readthedocs.io/en/latest/syntax/optional.html#math-shortcuts)
- [Linkify (with fuzzy matching disabled)](https://myst-parser.readthedocs.io/en/latest/syntax/optional.html#linkify)
- [Substitutions](https://myst-parser.readthedocs.io/en/latest/syntax/optional.html#substitutions-with-jinja2)
- [Auto-Generated Header Anchors](https://myst-parser.readthedocs.io/en/latest/syntax/optional.html#auto-generated-header-anchors)
- [Definition Lists](https://myst-parser.readthedocs.io/en/latest/syntax/optional.html#definition-lists)
- [Task Lists](https://myst-parser.readthedocs.io/en/latest/syntax/optional.html#task-lists)
- [Field Lists](https://myst-parser.readthedocs.io/en/latest/syntax/optional.html#field-lists)
- [Attributes](https://myst-parser.readthedocs.io/en/latest/syntax/optional.html#attributes)

Markdown documents may start with an optional YAML metadata header with any of the following configured.

```yaml
---
tocdepth: 2 # maximum depth for the page table of contents on the right sidebar
orphan: true # must be set if page is not included in the structure configuration
no-search: true # omit the page from text-based search
html_theme.sidebar_secondary.remove: true # remove the right sidebar for the page
---
```

Although possible, the inclusion of raw HTML within Markdown documents is strongly discouraged and should only be done to implement advanced functionality or accessibility improvements that otherwise would not be possible.

<!-- prettier-ignore -->

> [!NOTE]
> README documents should be written using [GitHub Flavored Markdown (GFM)](https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax) instead of MyST Markdown.

#### reStructuredText

reStructuredText (RST) is a plaintext markup language used primarily for technical documentation within the Python community. It is the default markup language in Sphinx and hence also supported by this project. However, the use of reStructuredText is discouraged and any new material should be written using MyST Markdown whenever possible. RST support exists primarily to allow the inclusion of preexisting legacy materials. Here are some resources to assist in understanding RST syntax and converting existing RST documents to MyST Markdown.

- [Sphinx reStructuredText Guide](https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html)
- [Migrating from reStructuredText to MyST Markdown](https://docs.readthedocs.com/platform/stable/guides/migrate-rest-myst.html)

### Executable

Executable content can be written using a variety of notebook formats described below. Executable code snippets or cells included in the notebooks are executed during the build process and any outputs are automatically included into the built document. Functionality to make the code snippets interactive and executable by the end user within the browser might be included in the future. All executable content regardless of format must adhere to the following guidelines.

- **All dependencies must be included in the `guides` environment.** Code included in the notebooks may not attempt to install any packages and all dependencies must be listed in `environment.yaml` with pinned versions. Ensure the functionality of the environment and the success of the build process after adding or updating any dependencies. Ensure added dependencies are easily removable by adding them in a single descriptive commit or using comments within `environment.yaml` to denote which dependencies are specific to a given notebook.
- **All data should be pre-downloaded and included within the `source` directory.** Unless specifically demonstrating the downloading of data, any data used within the notebooks should be pre-downloaded and included either within the same directory as the notebook or within a designated `data` directory that is a child of the notebook directory. Direct download links to any included data files should be included in the notebook to allow the end user to manually download and explore the data. Data files should not exceed 50 megabytes.
- **Computation should be relatively fast and lightweight.** Any notebook run should not exceed two minutes and use no more than two gigabytes of RAM. Note that these are upper limits -- faster and less-intensive computation is strongly encouraged. Remember that all notebooks get executed during every build and the functionality to run the code within the browser is being added. Hence this platform should be used for lightweight examples and demonstrations requiring significant computation should be hosted elsewhere.

#### Jupyter Notebook

Jupyter Notebooks using either the interactive Python kernel (`python3`) or the interactive R kernel (`ir`) are supported. [Cell magic commands](https://ipython.readthedocs.io/en/stable/interactive/magics.html) are allowed but discouraged as any OS-dependent functionality is not guaranteed. Markdown cells can include MyST Markdown syntax and the JupyterLab installation included in the project environment is equipped with an extension that adds MyST Markdown rendering support. Hence it is recommended to use the JupyterLab installation included in the project environment to develop notebooks. JupyterLab can be launched as follows.

```
jupyter lab
```

#### MyST Markdown Notebook

MyST Markdown includes the functionality for text-based Jupyter notebooks via the [MyST NB](https://myst-nb.readthedocs.io/en/latest/index.html) extension. These are written entirely in Markdown and include [special YAML metadata](https://myst-nb.readthedocs.io/en/latest/authoring/text-notebooks.html#notebook-level-metadata) that define the computation kernel and [special code cell directives](https://myst-nb.readthedocs.io/en/latest/authoring/text-notebooks.html#syntax-for-code-cells) that specify which code should be executed during runtime. MyST Markdown notebooks are very similar to R Markdown notebooks but support all MyST Markdown syntax. Both the interactive Python kernel (`python3`) and the interactive R kernel (`ir`) can be used. See below for resources.

- [MyST NB Text-Based Notebooks Guide](https://myst-nb.readthedocs.io/en/latest/authoring/text-notebooks.html)

#### R Markdown Notebook

R Markdown notebooks are supported with certain limitations. Package installations via CRAN are not allowed (all dependencies must be installed via `environment.yaml`) and Python is not supported. Syntax specific to R Markdown can be used but note that some functionality might be lost during the build process. R Markdown is not directly supported and thus all R Markdown documents are converted to Jupyter Notebooks using [JupyText](https://jupytext.readthedocs.io/en/latest/) and then executed using the interactive R kernel (`ir`). Hence certain elements might not be rendered as expected. Also note that any functionality specific to RStudio and the `knitr` package is not supported as the R Markdown documents are never actually _knit_. See below for resources.

- [R Markdown Guide](https://bookdown.org/yihui/rmarkdown/)
- [Example of a Converted R Markdown File](https://myst-nb.readthedocs.io/en/latest/authoring/custom-formats.html)

<!-- prettier-ignore -->

> [!CAUTION]
> The properly capitalized `Rmd` extension must be used for R Markdown files to be correctly identified.

## Style Guidelines

### Prose

Style guidelines for prose are in active development and subject to change. The following are current recommendations that are not enforced.

- Use title case for all headings.
- Do not use "you" or "we" and avoid addressing the reader directly.
- Use active voice for specific instructions and passive voice otherwise.
- Define any acronyms the first time they are mentioned.
- Numbers up to ten should be written using words instead of numerals.
- Avoid duplication of effort -- link out to existing internal or external materials whenever possible.
- Take care when linking out to external materials -- prefer official or non-commercial resources whenever possible.
- Use [substitutions](https://myst-parser.readthedocs.io/en/latest/syntax/optional.html#substitutions-with-jinja2) for frequently included text subject to change like emails or URLs.

<!-- prettier-ignore -->

> [!WARNING]
> **Do not include any copyrighted material.** Material under copyright protections intended to be used for educational content under _fair use_ and the _TEACH Act_ should be limited to Tufts affiliates only and not be included in any publicly accessible resource.

### Code

Python code should conform to the [Black](https://github.com/psf/black) style and R code should conform to the [Tidyverse](https://style.tidyverse.org/) style. Style guidelines are enforced and automatically applied via pre-commit hooks during commits (if configured) and pull requests. Manual effort to conform to style guidelines is not needed. Style guides for other languages might be added in the future.

## Accessibility

All content should be mindful of screen readers and color contrast guidelines. All images should have alternative text and any decorative elements be marked accordingly to be ignored by screen readers. The colors used by the [PyData Sphinx Theme](https://pydata-sphinx-theme.readthedocs.io/en/latest/user_guide/accessibility.html) have been carefully selected to meet accessibly guidelines and hence should not be manually modified. Further accessibility guidelines are in active development and to be added in the future.

## Subject Tags

Tags can be defined using the `tags` field in the file-wide metadata. The field content must be a single string representing a space-delimited list of tags. Tags can only contain lowercase letters, numbers, and hyphens (`-`) with no more than one consecutive hyphen. This is enforced and improperly formatted tags will result in an extension error during the build process.

### MyST and R Markdown

Tags can be defined in the YAML metadata header of the file as follows.

```yaml
---
tags: tag tag2 another-tag
---
```

### reStructuredText

Tags can be specified in the metadata field list at the top of the file as follows.

```rst
:tags: tag tag2 another-tag
```

### Jupyter Notebook

Tags can be added to the notebook metadata JSON as follows. The metadata JSON can be accessed via the Property Inspector in the top-right of the JupyterLab interface (gear icon) or by opening the notebook as a text document and locating the `"metadata"` field (usually located after the `"cells"` field) in the notebook JSON.

```json
{
  "tags": "tag tag2 another-tag"
}
```

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

### Updating the Environment

If there have been any updates to the `environment.yaml` environment specification file, the environment should be updated as follows.

```
mamba env update --file environment.yaml
```

In some cases a full rebuild of the environment is needed. The environment can be removed and recreated as follows.

```
mamba env remove --name guides --yes
mamba env create --file environment.yaml --yes
```

### Committing Changes

Frequent commits are encouraged. The subject line of the commit message should be informative and not exceed 50 characters. Commit message bodies are not required but encouraged to provide additional detail where needed. Avoid using sentence-style capitalization and punctuation in the commit message subject line. When using pre-commit hooks, avoid passing the commit message directly via the `-m` flag in case checks fail and automatic fixes are applied.

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

A manual review is not required when merging to the `develop` branch and the pull request can be merged when all automated status checks pass. Run the following to have this automatically happen. Manual intervention is needed when status checks fail and automatic fixes are not possible. Ask for help if needed.

```
gh pr merge --auto
```

The merging of the PR triggers a build and deployment of the development version of the website. Note that this could take several minutes. Once the updated website has been deployed, make sure to navigate to the [development version of the website][guides-dev-url] and ensure all changes are reflected as expected.

## Publishing Workflow

The following workflow incorporates updates from the `develop` branch into the published website. This is automatically triggered every week (if updates are detected) but can be manually triggered as follows anytime when updates to the published site are needed.

### Preparing for Publishing

Switch to the `staging` branch to prepare for publishing.

```
git switch staging
git pull
```

Merge the latest version of the `develop` branch into `staging` to incorporate the changes.

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

Merges to `main` require approval and a manual review. Adding a comment tagging a reviewer to grab their attention is encouraged. Comments can be added via the GitHub CLI as follows.

```
gh pr comment
```

Comments can also be added using the web interface which can be easily accessed as follows.

```
gh pr view --web
```

### Merging the Pull Request

The PR can be merged once it is manually approved and all checks pass. Run the following to have this automatically happen. Merging of the PR will trigger the publishing workflow and result in a new build of the public website.

```
gh pr merge --auto
```

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
