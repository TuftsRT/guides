########
 Fasta3
########

**************
 Introduction
**************

Fasta3 is a suite of programs for searching nucleotide or protein
databases with a query sequence.

**********
 Versions
**********

-  36.3.8

**********
 Commands
**********

-  fasta36
-  fastf36
-  fastm36
-  fasts36
-  fastx36
-  fasty36
-  ggsearch36
-  glsearch36
-  lalign36
-  ssearch36
-  tfastf36
-  tfastm36
-  tfasts36
-  tfastx36
-  tfasty36

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
   #SBATCH --job-name=fasta3      # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load fasta3/XXXX        ### Latest version is recommended.
