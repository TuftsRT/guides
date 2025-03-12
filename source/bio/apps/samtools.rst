##########
 Samtools
##########

**************
 Introduction
**************

Samtools is a set of utilities for the Sequence Alignment/Map (SAM)
format.

**********
 Versions
**********

-  1.16
-  1.17

If you require a version newer than these two, please visit the
`Samtools download and installation page
<https://www.htslib.org/download/>`_

As of 2024-12-10, the latest version available is 1.21.

**********
 Commands
**********

-  ace2sam
-  blast2sam.pl
-  bowtie2sam.pl
-  export2sam.pl
-  fasta-sanitize.pl
-  interpolate_sam.pl
-  maq2sam-long
-  maq2sam-short
-  md5fa
-  md5sum-lite
-  novo2sam.pl
-  plot-ampliconstats
-  plot-bamstats
-  psl2sam.pl
-  sam2vcf.pl
-  samtools
-  samtools.pl
-  seq_cache_populate.pl
-  soap2sam.pl
-  wgsim
-  wgsim_eval.pl
-  zoom2sam.pl

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
   #SBATCH --job-name=samtools    # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge
   module load samtools/XXXX ### you can run *module avail samtools* to check all available versions

**********************************
 Using Multi-Threads for Samtools
**********************************

Samtools supports multi-threading to speed up operations like sorting,
indexing, and depth calculation. Below are examples of common tasks:

Setting Up SLURM for Multi-Threading:

.. code:: bash

   #!/bin/bash
   #SBATCH --job-name=samtools_job
   #SBATCH --output=samtools.%j.out
   #SBATCH --error=samtools.%j.err
   #SBATCH --ntasks=1               # Number of tasks (keep as 1 for multi-threaded Samtools)
   #SBATCH --cpus-per-task=8        # Number of CPU cores per task
   #SBATCH --mem=16G                # Memory per node
   #SBATCH --time=04:00:00          # Time limit (hh:mm:ss)
   #SBATCH -p batch

   # Load Samtools module
   module load samtools/1.17

   # Navigate to the working directory
   cd $SLURM_SUBMIT_DIR

   # Run Samtools command
   samtools sort -@ $SLURM_CPUS_PER_TASK -o sorted.bam input.bam

*****************
 Reference links
*****************

`Samtools documentation <https://www.htslib.org/doc/samtools.html>`_

`Samtools tutorials from other resources
<https://davetang.org/wiki/tiki-index.php?page=SAMTools>`_
