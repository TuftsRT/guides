############
 Parabricks
############

**************
 Introduction
**************

NVIDIA's Clara Parabricks brings next generation sequencing to GPUs,
accelerating an array of gold-standard tooling such as BWA-MEM, GATK4,
Google's DeepVariant, and many more. Users can achieve a 30-60x
acceleration and 99.99% accuracy for variant calling when comparing
against CPU-only BWA-GATK4 pipelines, meaning a single server can
process up to 60 whole genomes per day. These tools can be easily
integrated into current pipelines with drop-in replacement commands to
quickly bring speed and data-center scale to a range of applications
including germline, somatic and RNA workflows.

**********
 Versions
**********

-  4.0.0-1
-  4.2.1-1

**********
 Commands
**********

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
   #SBATCH --job-name=parabricks  # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load parabricks/XXXX    ### Latest version is recommended.
