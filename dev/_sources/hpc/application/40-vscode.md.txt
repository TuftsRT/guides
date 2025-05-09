# VSCode

In this walkthrough, we'll examine how to set up a tunnel between the Cluster and your local installation of VSCode. This can be useful when you need to create, edit and publish sophisticated codebases on the Cluster. We can also use this method for running `Jupyter Notebooks` from the Cluster.

This process requires that you already have a Github account. Please feel free to make one [here](https://github.com/).

## Local Extensions (Optional)

These extensions need to be installed and system meet all prerequisites indicated on the official extension page. These are optional, but they add additional functionality which may be important for data science workflows.

### Example extensions

Here are some examples of useful local extensions:

[Remote-SSH](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-ssh)

[Remote Explorer](https://marketplace.visualstudio.com/items?itemName=ms-vscode.remote-explorer)

[Remote-Tunnels](https://marketplace.visualstudio.com/items?itemName=ms-vscode.remote-server)

[Jupyter](https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter) (optional)

## SSH Keyless Access (Optional)

Please follow the instructions here to setup your SSH keyless access. This will make your life much easier in the long run: [SSH Keyless Access](https://www.tecmint.com/ssh-passwordless-login-using-ssh-keygen-in-5-easy-steps/)

Once you are ready, let's get started.

## Tunnels:

- [Tmux](30-tmux) session (Optional)

Start a tmux session on Tufts HPC cluster in any shell environment on the login node.

- Allocate resources on HPC Cluster

Allocate appropriate amount of resources you need for your session with `srun` to start an [interactive session](../slurm/interactive.md) inside the tmux session.

> e.g. `srun -p interactive -n 2 --mem=4g --pty bash`

- Load vscode_cli module

`module load vscode_cli/1.77.3`

Then configure and start tunnel:

`code tunnel`

You will need to follow any Two Factor Authentication steps from Github to proceed. Once you have done so, copy the link given into your local browser. You should now see a VSCode window running from the browser. Feel free to connect any directory by clicking on the file explorer on the left. Currently, VSCode does not support Python environments to be ported through the remote tunnel. Read more [here](https://github.com/microsoft/vscode-python/issues/21557).

You can either use the VSCode server browser tab from there, or you can go back to your locally installed VSCode. You can find your active tunnels in "Remote Explorer".
