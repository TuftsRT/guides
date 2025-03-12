#########
 Bowtie2
#########

**************
 Introduction
**************

Bowtie 2 is an ultrafast and memory-efficient tool for aligning
sequencing reads to long reference sequences. It is particularly good at
aligning reads of about 50 up to 100s or 1,000s of characters, and
particularly good at aligning to relatively long (e.g. mammalian)
genomes. Bowtie 2 indexes the genome with an FM Index to keep its memory
footprint small: for the human genome, its memory footprint is typically
around 3.2 GB. Bowtie 2 supports gapped, local, and paired-end alignment
modes.

**********
 Versions
**********

-  2.4.2
-  2.5.1

**********
 Commands
**********

-  bowtie2
-  bowtie2-build
-  bowtie2-inspect

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
   #SBATCH --job-name=bowtie2     # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load bowtie2/XXXX       ### Latest version is recommended.
