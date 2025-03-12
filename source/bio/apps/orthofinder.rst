#############
 Orthofinder
#############

**************
 Introduction
**************

Orthofinder is a fast, accurate and comprehensive platform for
comparative genomics. It finds orthogroups and orthologs, infers rooted
gene trees for all orthogroups and identifies all of the gene
duplication events in those gene trees. It also infers a rooted species
tree for the species being analysed and maps the gene duplication events
from the gene trees to branches in the species tree. OrthoFinder also
provides comprehensive statistics for comparative genomic analyses.

**********
 Versions
**********

-  2.5.5

**********
 Commands
**********

-  orthofinder

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
   #SBATCH --job-name=orthofinder # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load orthofinder/XXXX   ### Latest version is recommended.
