##################
 Nf-core-ampliseq
##################

**************
 Introduction
**************

nfcore/ampliseq is a bioinformatics analysis pipeline used for amplicon
sequencing, supporting denoising of any amplicon and supports a variety
of taxonomic databases for taxonomic assignment including 16S, ITS, CO1
and 18S. Phylogenetic placement is also possible. Multiple region
analysis such as 5R is implemented. Supported is paired-end Illumina or
single-end Illumina, PacBio and IonTorrent data. Default is the analysis
of 16S rRNA gene amplicons sequenced paired-end with Illumina.

**********
 Versions
**********

-  2.8.0
-  2.9.0
-  2.10.0
-  2.11.0

**********
 Commands
**********

-  ampliseq

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
   #SBATCH --job-name=nf-core-ampliseq    # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load nf-core-ampliseq/XXXX      ### Latest version is recommended.
