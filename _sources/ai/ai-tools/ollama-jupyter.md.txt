# Using the Ollama + Jupyter Application on the Cluster

By Peter Nadel, Digital Humanities Natural Language Processing Specialist

In this document, we show the best practices for accessing Ollama with Jupyter on the Tufts High Performance Compute (HPC) Cluster. This document does not give an in-depth description of all of the features of Ollama or Jupyter. For this, we encourage you to explore the documentation of each of these projects:

- [Ollama documentation](https://docs.ollama.com/)
- [Jupyter documentation](https://jupyter.org/).

We also encourage you to explore the [HPC Cluster documentation](https://rtguides.it.tufts.edu/hpc/index.html).

## What is Ollama

Ollama is a program for using large language models (LLMs). Specifically, Ollama focuses open source and open weight LLMs usually accessible on websites like [HuggingFace](huggingface.co). We take advantage of the compute resources on our Cluster to run Ollama with a variety of pre-downloaded models. Ollama is the bridge that we use to send your questions/requests/responses from your keyboard to the LLM itself.

## What is Jupyter

Jupyter is a web-based Python programming platform that enables interactive computing through computational notebooks. A notebook is a shareable document that combines computer code, natural language documentation, interactive visualizations and data. Jupyter notebooks are designed for fast prototyping and code explanation. It is an ideal choice for classroom learning, as well as research software development.

## Who is this for?

This application is best suited for individuals with experience with the Python programming language and who would like to learn more about how AI can facilitate certain research applications. This tools can support high-throughput methods, but uses may find other approaches easier when scaling to larger and larger datasets. For more information on other options you have for AI on the Cluster, please visit [this link](TBD).

## Getting started

Ollama + Jupyter is an Open OnDemand application on the Cluster, meaning that it can be accessed from the Interactive Apps drop down in the Open OnDemand website. To get started, visit and log into the [Open OnDemand website for the Tufts Cluster](ondemand-p01.pax.tufts.edu). Once there, select the "Interactive Apps" drop down and click on "Ollama + Jupyter".

![Interactive Apps dropdown](../../images/OllamaJupyterList.png)

## Configuring your session

Once you've clicked on "Ollama + Jupyter", you will be able to configure the setting for using the application. Some of these options can be confusing, so we have left an example configuration below. If you are unsure, feel free to use this one. Otherwise, we explore what these parameters mean here:

- *Number of hours*: This parameter controls how long your session will run for. At the conclusion of this time, your session will end. Be sure to choose a time that matches how you expect to need in hours. You can always budget more time than you may need and they end the session early if you need.
- *Number of cores*: This field controls how many CPU cores are allocated for your session. It is important to pick a value proportional to the size of the LLM you'd like to run. If you are having trouble choosing, you can use the value shown below.
- *Amount of Memory (GB)*: This setting controls how many gigabytes of RAM are allocated to your session. This value can also be difficult to choose, so I like to use double the amount of CPU cores that I have selected.
- *Partition*: For most Ollama application, you should choose the "gpu" option. Generally, we require hardware acceleration to run LLMs. You can run some models, however, with just CPUs, especially if you adjust the number of cores and amount of memory to be quite high, in which case, you could select "batch" for this option.
- *GPU architecture*: This parameter controls the type of GPU that is allocated for your session. For the most part, it may not matter, however, if you pick a GPU type that is high demand, it may take longer for your session to get allocated. For more information on how to check demand for GPUs, use the `hpctools` CLI. Learn more [here](https://rtguides.it.tufts.edu/hpc/examples/hpctools.html).

The rest of the the fields should remain in their default configuration. When you are ready, click "Launch".

![Demo Config](../../images/OllamaJupyterConfig.png)

## Logging in to Ollama

Once you've launched your session, you will see the loading message below. It is very normal to see this for a couple minutes.

![Starting](../../images/OllamaJupyterStarting.png)

When the application is ready, you'll see the message below:

![Running](../../images/OllamaJupyterRunning.png)

When you are ready click "Connect to Jupyter". All of the data that you enter here **will stay on the Cluster and will not be shared on the internet.**

For a list of all models and research use cases can be found [here](<>). For a list of features, please visit this page: https://docs.jupyter.org/en/latest/. For other questions, please reach out to Research Technology at: tts-research@tufts.edu.
