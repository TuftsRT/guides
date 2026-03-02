# Batch Jobs

A Batch jobs is one where once submitted to the scheduler, the user does not intervene. The application is started by the scheduler once the requested nodes are allocated, and it runs unattended until completion. After completion the user returns at their convivence to review the program output files. Batch jobs are created by the user submitting a batch script file to the Slurm scheduler using the sbatch command.

- It is easy to submit a batch job.

- You will need a slurm batch job script which contains all the resource requirements and the commands to run your program.

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
#SBATCH --mail-user=Your_Tufts_Email@tufts.edu  #please make sure to use your own Tufts email

# [commands_you_would_like_to_exe_on_the_compute_nodes]
# have a clean start. purge all loaded modules in current environment
module purge
# for example, running a python script
# load the module so the correct version python is available to you
module load miniforge/25.3.0
# If you have a conda env that you would like to use, activate it 
source activate [target_env]
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
#SBATCH --mail-user=Your_Tufts_Email@tufts.edu #please make sure to use your own Tufts email

# [commands_you_would_like_to_exe_on_the_compute_nodes]
# have a clean start. purge all loaded modules in current environment
module purge
# for example, running a python script
# load the module so the correct version python is available to you
module load miniforge/25.3.0
# when using GPUs, make sure to load the appropriate version of cuda toolkit to provide necessary libraries for your application
module load cuda/12.9.0
# If you have a conda env that you would like to use, activate it 
source activate [target_env]
# run python script
python myscript.py #make sure myscript.py exists in the current directory
# make sure you save all plots, data, outputs generated to files in your script
# Don't forget to deactivate your conda env if you are using one
conda deactivate
```

**Submit** the job using the following command from command line interface:

`$ sbatch mygpujob.sh`

This will output the unique ID for the job. Keep this for reference.
