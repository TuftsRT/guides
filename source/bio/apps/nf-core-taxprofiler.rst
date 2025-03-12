#####################
 Nf-core-taxprofiler
#####################

**************
 Introduction
**************

nf-core/taxprofiler is a bioinformatics best-practice analysis pipeline
for taxonomic classification and profiling of shotgun short- and
long-read metagenomic data. It allows for in-parallel taxonomic
identification of reads or taxonomic abundance estimation with multiple
classification and profiling tools against multiple databases, and
produces standardised output tables for facilitating results comparison
between different tools and databases.

**********
 Versions
**********

-  1.1.5
-  1.1.6
-  1.1.7
-  1.1.8

**********
 Commands
**********

-  taxprofiler

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
   #SBATCH --job-name=nf-core-taxprofiler # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load nf-core-taxprofiler/XXXX   ### Latest version is recommended.
