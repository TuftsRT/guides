##########
 Fasttree
##########

**************
 Introduction
**************

FastTree infers approximately-maximum-likelihood phylogenetic trees from
alignments of nucleotide or protein sequences. FastTree can handle
alignments with up to a million of sequences in a reasonable amount of
time and memory. For large alignments, FastTree is 100-1,000 times
faster than PhyML 3.0 or RAxML 7.

**********
 Versions
**********

-  2.1.11

**********
 Commands
**********

-  fasttree
-  FastTree
-  FastTreeMP

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
   #SBATCH --job-name=fasttree    # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load fasttree/XXXX      ### Latest version is recommended.
