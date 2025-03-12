##########
 Bedtools
##########

**************
 Introduction
**************

Bedtools is an extensive suite of utilities for genome arithmetic and
comparing genomic features in BED format.

**********
 Versions
**********

-  2.30.0
-  2.31.0

**********
 Commands
**********

-  annotateBed
-  bamToBed
-  bamToFastq
-  bed12ToBed6
-  bedpeToBam
-  bedToBam
-  bedToIgv
-  bedtools
-  closestBed
-  clusterBed
-  complementBed
-  coverageBed
-  expandCols
-  fastaFromBed
-  flankBed
-  genomeCoverageBed
-  getOverlap
-  groupBy
-  intersectBed
-  linksBed
-  mapBed
-  maskFastaFromBed
-  mergeBed
-  multiBamCov
-  multiIntersectBed
-  nucBed
-  pairToBed
-  pairToPair
-  randomBed
-  shiftBed
-  shuffleBed
-  slopBed
-  sortBed
-  subtractBed
-  tagBam
-  unionBedGraphs
-  windowBed
-  windowMaker

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
   #SBATCH --job-name=bedtools    # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load bedtools/XXXX      ### Latest version is recommended.
