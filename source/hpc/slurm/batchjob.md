# Batch Jobs

Batch jobs are best suited for production level work, as there's no chance for users to interact and modify the programs while it's running. So be sure to debug and test your programs beforehand.

It is easy to submit a batch job.

You will need a slurm batch job script which contains all the resource requirements and the commands to run your program.

We strongly recommend user check available resources using [hpctools - Tufts HPC Helper Tool](../examples/hpctools.md) before submitting batch jobs.

Please see examples below:

## CPU Batch Job

Write a batch submission script e.g. **mycpujob.sh**

```bash
#!/bin/bash -l
#SBATCH -J My_Job_Name   #job name
#SBATCH --time=00-01:20:00  #requested time (DD-HH:MM:SS)
#SBATCH -p batch,preempt    #running on "batch" or "preempt" partition, wherever resource is available first
#SBATCH -N 1    #1 nodes #for many shared-memory programs,please leave -N as 1.
#SBATCH -n 2   #2 tasks total and 1 cpu per task, that gives you 2 cpu cores for this job
#SBATCH --mem=2g  #requesting 2GB of RAM total for the number of cpus you requested
#SBATCH --output=MyJob.%j.%N.out  #saving standard output to file, %j=JOBID, %N=NodeName
#SBATCH --error=MyJob.%j.%N.err   #saving standard error to file, %j=JOBID, %N=NodeName
#SBATCH --mail-type=ALL    #email options
#SBATCH --mail-user=Your_Tufts_Email@tufts.edu

#[commands_you_would_like_to_exe_on_the_compute_nodes]
# have a clean start. purge all loaded modules in current environment
module purge
# for example, running a python script
# load the module so the correct version python is available to you
module load anaconda/2021.05
# If you have a conda env that you would like to use, activate it here using "source activate xxx". DO NOT USE "conda activate"
source activate mycondaenv
# run python script
python myscript.py #make sure myscript.py exists in the current directory
# make sure you save all plots, data, outputs generated to files in your script
# Don't forget to deactivate your conda env if you are using one
conda deactivate
```

**Submit** the job using the following command from command line interface:

`$ sbatch mycpujob.sh`

You will receive an unique ID for the job.

## GPU Batch Job

Write a batch submission script e.g. **mygpujob.sh**

```bash
#!/bin/bash -l
#SBATCH -J My_Job_Name   #job name
#SBATCH --time=00-00:20:00  #requested time (DD-HH:MM:SS)
#SBATCH -p gpu,preempt    #running on "batch" or "preempt" partition, wherever resource is available first
#SBATCH -N 1    #1 nodes #for many shared-memory programs,please leave -N as 1.
#SBATCH -n 2   #2 tasks total and 1 cpu per task, that gives you 2 cpu cores for this job
#SBATCH --mem=2g  #requesting 2GB of RAM total for the number of cpus you requested
#SBATCH --gres=gpu:a100:1	#requesting 1 A100 GPU
#SBATCH --constraint="a100-80G" #only requesting the A100 GPU with 80GB memory
#SBATCH --output=MyJob.%j.%N.out  #saving standard output to file, %j=JOBID, %N=NodeName
#SBATCH --error=MyJob.%j.%N.err   #saving standard error to file, %j=JOBID, %N=NodeName
#SBATCH --mail-type=ALL    #email optitions
#SBATCH --mail-user=Your_Tufts_Email@tufts.edu

#[commands_you_would_like_to_exe_on_the_compute_nodes]
#
module purge
# for example, running a python script
# load the module so the correct version python is available to you
module load anaconda/2021.05
# when using GPUs, make sure to load the appropriate version of cuda toolkit to provide necessary libraries for your application
module load cuda/12.2
# If you have a conda env that you would like to use, activate it here using "source activate xxx". DO NOT USE "conda activate"
source activate [target_env]
# run python script
python myscript.py #make sure myscript.py exists in the current directory
# make sure you save all plots, data, outputs generated to files in your script
# Don't forget to deactivate your conda env if you are using one
conda deactivate
```

**Submit** the job using the following command from command line interface:

`$ sbatch mygpujob.sh`

You will receive an unique ID for the job.

\*\*You can find more sample scripts in `/cluster/tufts/hpc/tools/slurm_scripts`. \*\*

**Feel free to cope the scripts to your own directory and modify them for your own use.**
