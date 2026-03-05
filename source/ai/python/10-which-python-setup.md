# Which Python Setup is Right for You?

There are many ways to access or install Python, and many users (not just beginners!) may struggle to decide which setup to use and how to implement it. In this guide, we try to make it easier for you by presenting our recommendations for setting up Python if you're planning to use it specifically for statistics, research, data analysis, data visualization, or data science.

Note that this is far from an exhaustive list of possible setups. The options we've selected represent reliable, well-tested options for users with various experience levels and research needs. If you don't see your favorite IDE or package manager here[^1], don't worry -- you can go with what works for you. Just be sure to do the research on your setup for any possible licensing restrictions or security vulnerabilities.

## Our Recommendations

- **First-Time Python Learners and Beginner to Intermediate Users:**
  - [Google Colab](#google-colab-a-solid-option-for-beginner-to-intermediate-users)
- **Intermediate to Advanced Users, or Users with Sensitive/IRB-Restricted Data:**
  - [Microsoft Visual Studio Code with Conda (using the Miniforge distribution)](#microsoft-visual-studio-code-and-miniforge-conda-a-free-versatile-setup-for-intermediate-to-advanced-users)
- **Tufts High-Performance Computing (HPC) Users:**
  - Best for new HPC users and general interactive use: [Jupyter Lab Server](#best-for-new-hpc-users-and-general-interactive-use-jupyter-lab-server-app-on-the-hpc) Interactive App on the HPC
  - Recommendation for advanced users and batch jobs: [Microsoft VS Code Server](#recommendation-for-advanced-users-and-batch-jobs-microsoft-vs-code-server-interactive-app-and-conda-on-the-hpc) Interactive App on the HPC, in conjunction with conda on the HPC.
  - Best for purpose-built Python-based AI tools: see options [here](https://rtguides.it.tufts.edu/ai/ai-tools/index.html).

Below you'll find an overview of each of these options. When you have decided which setup is right for you, click on the links included in each section below to continue to the corresponding tutorials.

## Google Colab: A Solid Option for Beginner to Intermediate Users

Google Colab is a web-based service that allows you immediate access to Python without having to worry about the complexities of a local installation. It provides a user interface that is very similar to our recommended installation for intermediate/advanced users, making it easy to transfer your skills when you're ready to build a more advanced setup on your machine. Integration with Google Drive makes it easy to upload, organize, and analyze your data and files. The interface includes useful tools for coding such as syntax highlighting, linting (basic error checking), and code completion, and by default it also includes AI integration that can help you write your code and understand your errors. These features make Colab an ideal option for Python beginners.

> **Beginner Tip:** Many online courses and tutorials use Jupyter Notebook to teach Python. Because Colab is built on top of Jupyter Notebook, you can also use Colab for any tutorial or class that uses Jupyter Notebook to teach Python.

Google Colab is free to all users, but free users face restrictions on usage of computer resources. A premium version exists for users who want access to more memory and faster GPUs. Fortunately, the premium version is free for students and academics. For more details on the benefits of premium, see this link:
https://colab.research.google.com/signup

Google Colab is powerful enough that you may be able to use it for all your research, with access to high-performance GPUs and an extensive collection of pre-installed libraries. (You can also install additional libraries with pip as needed.) It also makes for easy collaboration, since everything is online and can be accessed from any web browser with no setup.

Colab has its limits, however. First and foremost, it should not be used with sensitive or IRB-restricted data (see the note below). Second, due to session time limits and restrictions on heavy use, it is not ideal for large projects that may require long computing times. And finally, because it does not provide a beginner-friendly way to manage Python environments, older code may stop working when it becomes incompatible with newer versions of Python and various libraries, making it hard to reuse code or revisit older projects. For these reasons, we only recommend Google Colab for beginners; intermediate to advanced users will want to proceed to the next section.

> **Do not use Google Colab with IRB-restricted or sensitive data.** To analyze data in Colab, data files must first be uploaded to Google Drive, which is not approved by many Data Use Agreements (DUAs) or IRB policies. Researchers working with such data should use a local installation (see the Microsoft Visual Studio Code and Miniforge section below) or consider applying for an account with the Tufts High-Performance Computing (HPC) cluster. If you're unsure about where to securely store your data, check the [data storage finder](https://access.tufts.edu/data-finder). For advanced consultations on data storage and computing needs, contact datalab-support@tufts.edu.

### Summary of Pros and Cons for Google Colab:

**Pros:**

- Zero setup: No local installs, instant Python environment
- Free to use: Students and academics can use either the free version or gain free access to the premium version
- Useful coding tools, including syntax highlighting, linting (basic error checking), code completion, and AI integration
- Extensive collection of pre-installed libraries, ability to install additional libraries with pip
- Google Drive integration: Simple file storage and sharing
- Cloud execution with access to powerful GPUs
- Interface comparable to more advanced installations: skills transfer easily, great for learners

**Cons:**

- Data privacy concerns: Often inappropriate for sensitive or IRB-restricted data
- Weak for long jobs: Not suitable for multi-day or large batch runs due to session time-outs and restrictions on heavy use
- No beginner-friendly way to manage Python environments: old code may become incompatible with newer versions of Python and its libraries
- (Advanced) Can only use PIP for installing packages, so may not be compatible with C libraries or certain command line tools that are needed for more advanced workflows.

If you think Google Colab is right for you, [click here for instructions on how to get started.](20-google-colab.md)

## Microsoft Visual Studio Code and Miniforge (Conda): A Free, Versatile Setup for Intermediate to Advanced Users

While Colab is a great option for Python learners and for many research purposes, intermediate to advanced users will often find the need to move to a local installation that allows them to carefully manage different versions of Python and its libraries, run code for longer sessions, and work with sensitive or IRB-protected data. (Note that even local storage may not be acceptable for some highly sensitive data; for more information, see the [Tufts' IRB documentation](https://viceprovost.tufts.edu/about-tufts-irb) or consult the [data storage finder](https://access.tufts.edu/data-finder).)

Our recommended local installation comes in two parts:

1. **Microsoft Visual Studio Code** is an IDE, or "Integrated Development Environment". An IDE provides a user interface where you can type your code, run it, and view its output.

1. **Miniforge** is a lightweight installer for **conda**, a tool for managing Python environments and libraries.

If you're new to programming, you might be surprised by the second part, and you may be wondering, _why do we need conda?_ If so, we recommend that you continue reading for an overview of Python environments and libraries. More advanced users may wish to skip ahead to our reasons for recommending Conda and VS Code, or to our summary of the pros and cons of this setup. Those who are ready may also go directly to our [installation and setup instructions.](30-vs-code-with-miniforge.md)

### What Are Environments and Libraries?

Python is constantly evolving, and new versions are released regularly. If you write code using one version of Python and later upgrade to a newer version, it's possible that your code may stop working or produce inconsistent results.

In addition, much of Python’s power comes from **libraries** (also called "packages"), which are bundles of third-party code written to perform specialized tasks like data analysis, visualization, or machine learning. These libraries are also continuously updated, and different versions of a library may depend on specific versions of Python or on other libraries (called "dependencies"). Mixing incompatible versions of Python, libraries, and dependencies can cause errors.

To manage this complexity, Python users create environments. An environment is an isolated workspace that contains:

- A specific version of Python
- A specific set of package versions

With an environment manager, individual projects can be given their own dedicated environments. This prevents updates in one project from breaking another and makes your work more stable and reproducible. An environment/package manager like **conda** helps you create these environments and install compatible versions of Python and its packages automatically.

### Why Conda?

**Conda** is an open-source tool for managing environments, Python versions, libraries, and dependencies. Though other tools are available, we recommend conda because it manages both Python itself and the libraries your projects depend on in a single system. Unlike tools that only install Python packages, conda also handles non-Python dependencies (such as compiled numerical libraries), which are common in data science and scientific computing. This makes installation more reliable, especially for packages like NumPy, SciPy, or PyTorch. Conda also makes it easy to create isolated environments with specific Python versions, reproduce those environments across machines, and avoid conflicts between projects. For research and data analysis workflows, this combination of flexibility, stability, and reproducibility makes conda a practical, powerful, and beginner-friendly choice.

Note that commonly available instructions for downloading and installing packages will often be based on the default pip command, as in "pip install pandas". This is because online documentation typically does not want to assume that you are using conda and want to make the instructions as general as possible. To gain full advantage of conda, however, you should get used to using conda for installing libraries instead. Thus a command for installing pandas would look like:

```
conda install pandas
```

For more information on conda vs pip, see https://www.anaconda.com/blog/understanding-conda-and-pip.

### What is Miniforge?

**Miniforge** is a lightweight installer for conda. While alternatives exist, we recommend Miniforge for its popularity, its balance of speed and ease of use, and for being fully open source and free for all users and use cases. Note that installing conda with Miniforge will require the user to manage environments and install packages via a terminal or command line. While some users may find this intimidating at first, there are relatively few commands that you have to learn, and the process quickly becomes easy with practice.

### Why Microsoft Visual Studio (VS) Code?

Visual Studio Code is a popular, cross-platform IDE (Integrated Development Environment) maintained and distributed by Microsoft. Though not fully open source, it is free to use for academic, personal, and commercial purposes.

VS Code works well with conda environments and provides helpful tools for Python programming, such as code completion, syntax highlighting and linting (basic error checking), integrated debugging with breakpoints, and built-in Git version control tools.

VS Code supports Jupyter Notebooks through the Jupyter extension and provides additional tools for viewing your data and the objects in your Jupyter environment. Alternatively, VS Code also allows you to use Jupyter-style interactive cells directly inside standard Python (.py) files using # %% cell markers. This is useful when you want the flexibility of interactive, cell-based execution during development, while still keeping your code in a clean, script-based format that is easier to version control, share, and run on remote systems.

More advanced users may also wish to use VS Code's Remote-SSH extension for access to the Tufts High-Performance Computing (HPC) cluster, which provides a useful way to edit files in the HPC filesystem and run your code with SLURM batch scripts.

### Summary of Pros and Cons for Microsoft VS Code with Miniforge (Conda)

**Pros**

- Free for all users and use cases
- Cross-platform (Windows, Mac OS, Linux)
- Explicit management of environments and libraries ensures stable and reproducible workflows
- Supports Jupyter Notebooks with additional tools for viewing your data and objects in your Jupyter environment; Jupyter cells can also be used within .py files
- Coding tools such as syntax highlighting, linting (basic error checking), and code completion; AI extensions available as well
- Integrated debugging with breakpoints
- Easy integration with Git version control
- Advanced: integrated remote SSH access to the Tufts HPC

**Cons**

- Initial setup and installation are more complex and can be confusing for beginners
- Initial learning curve: users must learn to manage environments and libraries through a terminal/command line using conda and access these environments from within VS Code

If you decide a local installation with Microsoft VS Code and Miniforge (Conda) is right for you, [click here for setup instructions.](30-vs-code-with-miniforge.md)

## Python on the Tufts High-Performance Computing Cluster (HPC)

If your research requires more computing power than your personal computer can provide, you may want to consider applying for an account to use the Tufts High-Performance Computing Cluster, or HPC.

For more information about using the HPC, see the HPC documentation [here](https://rtguides.it.tufts.edu/hpc/access/index.html).

There are several ways to program and run code on the HPC. We recommend the following:

### Best for new HPC users and general interactive use: Jupyter Lab Server app on the HPC

For those familiar with Jupyter Lab or Jupyter Notebooks, the Jupyter Lab Server app is the HPC implementation of that. Log in to the cluster [here](https://ondemand-p01.pax.tufts.edu/pun/sys/dashboard) and click on the "Interactive Apps" menu to see a list of all apps. Scroll down to the "Servers" section of the drop-down menu and select "Jupyter". You will be taken to a form that allows you to specify the resources you would like to request for your session. When you have specified your desired resources, click "Launch". Your session request will be placed in a queue until the resources are allocated, at which point a link will appear to access your session.

### Recommendation for advanced users and batch jobs: Microsoft VS Code Server Interactive App and conda on the HPC

The VS Code and conda set up that we recommend for intermediate/advanced local applications can be replicated on the HPC as well, with some additional steps and some slight differences in the user interface and workflow.

Note that the instructions for using this setup on the HPC are less detailed than the instructions we provide for the local installation version of this setup. If using conda and VS Code are unfamiliar to you and you want more detailed guidance on how to use them effectively, you may wish to follow our guides for [installing](30-vs-code-with-miniforge.md) and [using](40-VS-code-and-conda-users-guide.md) this setup on a local computer first.

#### Conda on the HPC

To learn how to manage conda environments on the HPC, click [here](https://rtguides.it.tufts.edu/hpc/application/10-condaenv.html).

#### VS Code Server

VS Code Server provides a user interface that is similar to what you would get from a local installation of VS Code, with some minor differences. More documentation on using VS Code Server will be coming soon, including information on how to connect VS Code to your HPC-based conda environments.

To access VS Code Server on the HPC, log in to the HPC [here](https://ondemand-p01.pax.tufts.edu/pun/sys/dashboard) and click on the "Interactive Apps" menu. Scroll down to the "Servers" section of the drop-down menu and select "VS Code Server". You will be taken to a form that allows you to specify the resources you would like to request for your session. When you have specified your desired resources, click "Launch". Your session request will be placed in a queue until the resources are allocated, at which point a link will appear to access your session.

### For Purpose-built Python-based AI tools:

The HPC has a number of purpose-built tools for AI tasks, including speech and text recognition, secure LLM access, etc. For more information on these options [here](https://rtguides.it.tufts.edu/ai/ai-tools/index.html).

[^1]: Some popular IDEs not listed here include Spyder, Jupyter Lab, and PyCharm. There are also many options for [data science notebooks](https://bpb-us-e1.wpmucdn.com/sites.tufts.edu/dist/8/3138/files/2022/10/DS-Notebook-Comparison-Data-Lab.pdf).
