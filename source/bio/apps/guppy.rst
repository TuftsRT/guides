#######
 Guppy
#######

**************
 Introduction
**************

Guppy is a data processing toolkit that contains the Oxford Nanopore
Technologies’ basecalling algorithms, and several bioinformatic
post-processing features.

**********
 Versions
**********

-  6.4.6
-  6.5.7

**********
 Commands
**********

-  guppy_aligner
-  guppy_barcoder
-  guppy_basecall_client
-  guppy_basecaller
-  guppy_basecaller_duplex
-  guppy_basecaller_supervisor
-  guppy_basecall_server

*************
 Example job
*************

Adjust slurm options based on job requirements (`slurm cheat sheet
<https://slurm.schedmd.com/pdfs/summary.pdf>`_):

.. code::

   #!/bin/bash
   #SBATCH -p partitionName  # batch, gpu, preempt, mpi or your group's own partition
   #SBATCH -t 1:00:00  # Runtime limit (D-HH:MM:SS)
   #SBATCH -N 1   # Number of nodes
   #SBATCH -n 1   # Number of tasks per node
   #SBATCH -c 4   # Number of CPU cores per task
   #SBATCH --mem=8G       # Memory required per node
   #SBATCH --job-name=guppy       # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load guppy/XXXX ### Latest version is recommended.
