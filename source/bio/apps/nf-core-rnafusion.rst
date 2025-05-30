###################
 Nf-core-rnafusion
###################

**************
 Introduction
**************

nf-core/rnafusion is a bioinformatics best-practice analysis pipeline
for RNA sequencing consisting of several tools designed for detecting
and visualizing fusion genes. Results from up to 5 fusion callers tools
are created, and are also aggregated, most notably in a pdf visualiation
document, a vcf data collection file, and html and tsv reports.

**********
 Versions
**********

-  3.0.1
-  3.0.2

**********
 Commands
**********

-  rnafusion

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
   #SBATCH --job-name=nf-core-rnafusion   # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load nf-core-rnafusion/XXXX     ### Latest version is recommended.
