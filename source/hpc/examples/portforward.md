# Port Forwarding on HPC Cluster

> In this walkthrough, we'll see how to run interactive applications from the Cluster on GPUs.
> This is particularly useful for running `Jupyter Lab` and `Jupyter Notebook`, but can will work for any applications that run from a port on a Cluster compute node.

## 1. Login to the Cluster from terminal with SSH

```
$ ssh your_username@login.pax.tufts.edu
```

Follow any Two-Factor Authentication directions that are necessary.

## 2. (Optional) Start a tmux session on the login node

Tmux allows us to run multiple terminals from the same interface. It makes the port forwarding process a bit simpler, but not required.

```
$ module load tmux
$ tmux new -s mynewsession
```

Further details on [Tmux](../application/30-tmux)

## 3. Start interactive session

More [GPU resource](../compute/gpu) details

First, we'll need to request GPU resources and then load in the required software for working with CUDA and Python.

```
$ srun -p gpu --gres=gpu:a100:1 -n 8 --mem=16g -t 3-0 --pty bash
$ module load cuda/12.2
$ module load anaconda/2021.05
```

More instructions on requesting resources can be find in the [Slurm](../slurm/index) section

You will need to activate the conda env you are interested in using and create a ipykernel for it. To continue, enter the following:

```
$ source activate
```

To check possible conda environments, use the following command:

```
$ conda env list
```

As an example, we'll use the `cupy` environment.

```
$ conda activate cupy
```

You should now see the name of your environment in parentheses in front of your command prompt:

```
$ (cupy) [your_username@your_compute_node ~]$
```

### Before we can run `Jupyter`, make sure to be in the directory you need to be and know the node name of your compute node

To isolate the name of your compute node, enter the following

```
$ cd ~
$ hostname | cut -f 1 -d .
```

### 4. Start `Jupyter Notebook` or `Jupyter Lab`

Make sure that you start it without a browser and you pick a port that is not used. You can expect to see similar output:

**For Jupyter Notebook:**

```
$ jupyter notebook --no-browser --port=7777

[I 15:41:01.399 NotebookApp] Serving notebooks from local directory: /cluster/home/your_username
[I 15:41:01.399 NotebookApp] Jupyter Notebook 6.4.12 is running at:
[I 15:41:01.399 NotebookApp] http://localhost:7777/?token=082b7ce58015a6189458316f4b75aa22d1b8b2108705afb1
[I 15:41:01.399 NotebookApp]  or http://127.0.0.1:7777/?token=082b7ce58015a6189458316f4b75aa22d1b8b2108705afb1
[I 15:41:01.399 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 15:41:01.404 NotebookApp]

    To access the notebook, open this file in a browser:
        file:///cluster/home/your_username/.local/share/jupyter/runtime/nbserver-49351-open.html
    Or copy and paste one of these URLs:
        http://localhost:7777/?token=082b7ce58015a6189458316f4b75aa22d1b8b2108705afb1
     or http://127.0.0.1:7777/?token=082b7ce58015a6189458316f4b75aa22d1b8b2108705afb1
```

**For Jupyter Lab:**

```
$ jupyter lab --no-browser --port=7777

[I 2022-10-05 15:53:26.897 ServerApp] jupyterlab | extension was successfully linked.
[I 2022-10-05 15:53:27.433 ServerApp] nbclassic | extension was successfully linked.
[I 2022-10-05 15:53:28.433 LabApp] JupyterLab extension loaded from /cluster/tufts/hpc/tools/anaconda/202105/lib/python3.8/site-packages/jupyterlab
[I 2022-10-05 15:53:28.433 LabApp] JupyterLab application directory is /cluster/tufts/hpc/tools/anaconda/202105/share/jupyter/lab
[I 2022-10-05 15:53:28.437 ServerApp] jupyterlab | extension was successfully loaded.
[I 2022-10-05 15:53:28.467 ServerApp] nbclassic | extension was successfully loaded.
[I 2022-10-05 15:53:28.468 ServerApp] Serving notebooks from local directory: /cluster/home/your_username
[I 2022-10-05 15:53:28.468 ServerApp] Jupyter Server 1.4.1 is running at:
[I 2022-10-05 15:53:28.468 ServerApp] http://localhost:7777/lab?token=07509cd7fb3e7dfcbee37fdd9514b96fe627e78d3178d4aa
[I 2022-10-05 15:53:28.468 ServerApp]  or http://127.0.0.1:7777/lab?token=07509cd7fb3e7dfcbee37fdd9514b96fe627e78d3178d4aa
[I 2022-10-05 15:53:28.468 ServerApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 2022-10-05 15:53:28.477 ServerApp]

    To access the server, open this file in a browser:
        file:///cluster/home/your_username/.local/share/jupyter/runtime/jpserver-50343-open.html
    Or copy and paste one of these URLs:
        http://localhost:7777/lab?token=07509cd7fb3e7dfcbee37fdd9514b96fe627e78d3178d4aa
     or http://127.0.0.1:7777/lab?token=07509cd7fb3e7dfcbee37fdd9514b96fe627e78d3178d4aa
```

## 5. Port Forwarding from the terminal

- **If you are using Tmux:**

1. Press `ctrl` + `b`
1. Press `c`
   This will create and open a new terminal window.

- **If you are not using Tmux:**
  Open a NEW terminal.

**For both:**
SSH into your compute node using port forwarding. Follow the model below. In this case, the `Jupyter Notebook` was started on node "s1cmp002" and port is "7777", but this may change for your case.

```
$ ssh -tt your_username@login.cluster.tufts.edu -L 7777:localhost:7777 ssh s1cmp002 -L 7777:localhost:7777
```

## 6. Using the application

Now open a browser window, and copy & paste in the URL:

```
http://localhost:7777/?token=082b7ce58015a6189458316f4b75aa22d1b8b2108705afb1
```

## Enjoy!

If you need something to play with, there are example notebooks with notes in `/cluster/tufts/hpc/tools/GPU_with_Python`. Please copy the folder or the notebooks to your personal space (home directory /cluster/home/your_username or research storage space /cluster/tufts/your_lab_name).
