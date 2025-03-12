#########
 Subread
#########

**************
 Introduction
**************

Subread carries out high-performance read alignment, quantification and
mutation discovery. It is a general-purpose read aligner which can be
used to map both genomic DNA-seq reads and DNA-seq reads. It uses a new
mapping paradigm called seed-and-vote to achieve fast, accurate and
scalable read mapping. Subread automatically determines if a read should
be globally or locally aligned, therefore particularly powerful in
mapping RNA-seq reads. It supports INDEL detection and can map reads
with both fixed and variable lengths.

**********
 Versions
**********

-  1.6.4
-  2.0.1

**********
 Commands
**********

-  detectionCall
-  exactSNP
-  featureCounts
-  flattenGTF

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
   #SBATCH --job-name=subread     # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load subread/XXXX       ### Latest version is recommended.
