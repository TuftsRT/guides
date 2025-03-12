# Cluster Login

```{important}
   VPN - For Off-campus access from non-Tufts Networks please connect to [the Tufts VPN](https://access.tufts.edu/vpn).
```

```{important}
   Requires an Active Cluster Account, to request an account [click here](10-account-request).
```

There are two main ways to access the HPC cluster to run software:

1. Web browser access via the [OnDemand Portal](https://ondemand.pax.tufts.edu)
1. Command line access via Secure Shell (SSH)

## New Users

For new users who have not used HPC resources before, we recommend starting with the **OnDemand Portal** via your web browser. OnDemand is useful for running many graphical applications such as MATLAB, Jupyter Notebooks, RStudio, and uploading downloading small files (\<900MB) all within your web browser session.

Users who require command line access, with greater control and longer running times for their research tasks are advised to use **Secure Shell (SSH)**.

## 1) OnDemand Portal

To access the OnDemand Portal, open your web browser. Note that for the best user experience we strongly recommended Chrome or Firefox. You cannot use the default MacOS browser, Safari.

Navigate to the Tufts HPC [OnDemand](https://ondemand.pax.tufts.edu). If you are off-campus, [the Tufts VPN](https://access.tufts.edu/vpn) is required.

You will see a login screen as below, and you should use your Tufts username (lowercase) and password to login:

![OnDemand Login Screen](https://tufts.box.com/shared/static/kwskz4ybo21sq2agcq4e7khmknc99qa0.jpg)

You'll see the OnDemand logo, along with tabs including:

- Files
- Jobs
- Clusters
- Interactive Apps
- Bioinformatics Apps
- ...

You can now use the Tufts HPC cluster resources on OnDemand.
[Getting Started with OnDemand](15-ondemand)

## 2) SSH

The SSH protocol aka **Secure Shell** is a method for secure remote login from one computer to another via the command line. The command line, also known as the command prompt or the terminal, allows you to control your computer directly, or to connect to other computers (like the HPC) via SSH. On Windows, this is commonly done with Command Prompt or Windows Terminal. On MacOS, this is done via Terminal or other shells, and on Linux there are a variety of built-in Terminal programs such as Gnome Terminal, Console, and XTerm.

SSH is predominantly a text-based login method, but **Graphical User Interface (GUI)** forwarding can be enabled via the X Window System (X11) to access graphical programs running on a cluster node. The X11 window system allows us to forward graphics through a Secure Shell (SSH) session.

```{important}
   SSH login host is login.pax.tufts.edu (formerly login.cluster.tufts.edu)
```

- Command for Shell environment (default: bash):

  - `$ ssh your_utln@login.pax.tufts.edu`

- Explanation of the Command:

  - `ssh`: The command to initiate an SSH connection.The SSH protocol aka **Secure Shell** is a method for secure remote login from one computer to another.

  - `your_utln@login.pax.tufts.edu`: Replace `your_utln` with your actual Tufts username. This is the address of the remote HPC cluster.

  - You can add the argument `-X` to your ssh command to enable the X11 graphical forwarding in your SSH connection. X Window System must be installed on your computer. That looks like this:
    `ssh -X your_utln@login.pax.tufts.edu`

- What Happens After Running the Command for Login:

  - You will be prompted to enter your password.
  - You will need to complete the 2FA process if required (Tufts Network).

Now you are on a **Login Node** of the cluster (login-prod-0x) and in your **Home Directory** (~ or _/cluster/home/your_utln_).

`$ [your_utln@login-prod-02 ~]`

There are 3 login nodes: `login-prod-01`, `login-prod-02`, `login-prod-03`

**_Please DO NOT run any software on the login nodes._**

- Recommended SSH Client software (Optional)

  > - Cross-Platform: [SecureCRT](https://access.tufts.edu/securecrt)
  > - Native: Windows 11, Mac and Linux all include a command line SSH clients in the operating system (OS). Commands such as `scp` can be used for quick connections without installing additional software.
  > - Windows: [PuTTY](https://www.putty.org/)

- Setting up [SSH keyless access](https://www.tecmint.com/ssh-passwordless-login-using-ssh-keygen-in-5-easy-steps/)

  - Be sure your `~/.ssh` permission is correct! Otherwise, SSH won't work properly.

  - `. ssh` directory: 700 ( drwx------ )

  - public key ( `. pub` file): 644 ( -rw-r--r-- )

  - private key (` id_rsa` ): 600 ( -rw------- )

```{important}
   Please DO NOT run any software on the login nodes.
```
