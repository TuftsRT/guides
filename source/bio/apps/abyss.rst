#######
 Abyss
#######

**************
 Introduction
**************

ABySS is a de novo sequence assembler intended for short paired-end
reads and genomes of all sizes.

**********
 Versions
**********

-  2.3.7

**********
 Commands
**********

-  ABYSS
-  abyss-align
-  abyss-bloom
-  abyss-bloom-dbg
-  abyss-bowtie
-  abyss-bowtie2
-  abyss-bwa
-  abyss-bwamem
-  abyss-bwasw
-  abyss-db-txt
-  abyss-dida
-  abyss-fac
-  abyss-fatoagp
-  abyss-filtergraph
-  abyss-fixmate
-  abyss-fixmate-ssq
-  abyss-gapfill
-  abyss-gc
-  abyss-index
-  abyss-junction
-  abyss-kaligner
-  abyss-layout
-  abyss-longseqdist
-  abyss-map
-  abyss-map-ssq
-  abyss-mergepairs
-  abyss-overlap
-  ABYSS-P
-  abyss-paired-dbg
-  abyss-paired-dbg-mpi
-  abyss-pe
-  abyss-rresolver-short
-  abyss-samtoafg
-  abyss-scaffold
-  abyss-sealer
-  abyss-stack-size
-  abyss-tabtomd
-  abyss-todot
-  abyss-tofastq
-  AdjList
-  Consensus
-  DAssembler
-  DistanceEst
-  DistanceEst-ssq
-  KAligner
-  konnector
-  logcounter
-  MergeContigs
-  MergePaths
-  Overlap
-  ParseAligns
-  PathConsensus
-  PathOverlap
-  PopBubbles
-  SimpleGraph

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
   #SBATCH --job-name=abyss       # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load abyss/XXXX ### Latest version is recommended.
