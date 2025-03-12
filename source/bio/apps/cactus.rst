########
 Cactus
########

**************
 Introduction
**************

Cactus is a reference-free whole-genome multiple alignment program.

**********
 Versions
**********

-  2.8.1
-  2.8.1-gpu

**********
 Commands
**********

-  cactus
-  cactus2hal.py
-  cactus2hal-stitch.sh
-  cactus-align
-  cactus-align-batch
-  cactus_analyseAssembly
-  cactusAPITests
-  cactus_barTests
-  cactus_batch_mergeChunks
-  cactus-blast
-  cactus_chain
-  cactus_consolidated
-  cactus_covered_intervals
-  cactus_fasta_fragments.py
-  cactus_fasta_softmask_intervals.py
-  cactus_filterSmallFastaSequences.py
-  cactus-graphmap
-  cactus-graphmap-join
-  cactus-graphmap-split
-  cactus_halGeneratorTests
-  cactus_local_alignment.py
-  cactus_makeAlphaNumericHeaders.py
-  cactus-minigraph
-  cactus-prepare
-  cactus-prepare-toil
-  cactus-preprocess
-  cactus-refmap
-  cactus_softmask2hardmask

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
   #SBATCH --job-name=cactus      # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load cactus/XXXX        ### Latest version is recommended.
