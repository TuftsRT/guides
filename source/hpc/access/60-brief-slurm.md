# Running your software

There are two ways to run software on the HPC Cluster, OnDemand ( Simplest ) and via SSH ( Most flexible ). Both
utilized the Slurm job scheduler to allocate you resources (CPU, Memory, and GPU if needed) on the compute nodes.
We will cover both briefly on this page. Each is also covered extensively in their own sections [Slurm Job Scheduler](../slurm/index) and [OnDemand](../access/15-ondemand) .

## OnDemand (Simplest)

OnDemand is a web portal that provides easy access to the HPC environment using a web browser. It provides web
forms for commonly used applications.

1. Start by going to [**OnDemand**](https://ondemand.pax.tufts.edu), https://ondemand.pax.tufts.edu using **Chrome or FireFox**
1. Login using your Tufts TTS username (utln) and password
1. **Interactive Apps** and **Bioinformatics Apps** menus provide a listing of the software available.
1. Select the applications of interest and fill in the appropriate values for your job.
1. Click the **Launch** button to submit your job to the Slurm scheduler. Your job will start once resources are
   available, this is typically within a few minutes, but could be longer during periods of high use.

Under **Jobs** -> **Active Jobs** you can see your running jobs. It also provides a file explorer Under **Files**.

For additional instructions on using it visit [OnDemand](../access/15-ondemand).

## SSH ( Most flexible )

Logging into a cluster login node using SSH and running a Slurm job is the most flexible way to utilize the
cluster. This method allows you to run any compatible software on the cluster. This includes custom software that
you wrote, or third party software you installed.

1. Login to login.pax.tufts.edu using SSH

1. Create a batch script file to run using Slurm , for example `nano myFirstJob.sh`

   ```
   #!/bin/bash -l
   #SBATCH -J My_First_Job     #job name
   #SBATCH --time=00-00:20:00  #requested time (DD-HH:MM:SS)
   #SBATCH -p batch
   #SBATCH -N 1    #1 nodes
   #SBATCH --output=MyFirstJob.%j.%N.out
   #SBATCH --error=MyFirstJob.%j.%N.err
   #SBATCH --mail-type=ALL    #email options
   #SBATCH --mail-user=Your_Tufts_Email@tufts.edu

   echo "My First Job"
   date

   ls #Put your program to run here

   echo "Done"
   ```

1. Submit your job using `sbatch myFirstJob.sh`

1. View your job in the queue using `squeue --Your_Tufts_Username`

For an in-depth description of all the Slurm options visit [Slurm Job Scheduler](../slurm/index).
