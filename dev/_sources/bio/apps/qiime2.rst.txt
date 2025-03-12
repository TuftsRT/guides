########
 Qiime2
########

**************
 Introduction
**************

QIIME 2 is a powerful, extensible, and decentralized microbiome analysis
package with a focus on data and analysis transparency. QIIME 2 enables
researchers to start an analysis with raw DNA sequence data and finish
with publication-quality figures and statistical results.

**********
 Versions
**********

-  2023.2
-  2023.5
-  2023.7
-  2023.9
-  2024.2

**********
 Commands
**********

-  biom
-  python
-  python3
-  qiime

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
   #SBATCH --job-name=qiime2      # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load qiime2/XXXX        ### Latest version is recommended.
