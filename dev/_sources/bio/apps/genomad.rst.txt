#########
 geNomad
#########

**************
 Introduction
**************

geNomad: Identification of mobile genetic elements.

Homepage: https://github.com/apcamargo/genomad

**********
 Versions
**********

-  1.8.1

**********
 Commands
**********

-  genomad

**********
 Database
**********

The genomad database has been downloaded for users. The current version
is **v1.7**. The database is stored in
`/cluster/tufts/biocontainers/datasets/genomad_db/v1.7`. Users can
follow the below example to run genomad with the downloaded database.

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
   #SBATCH --job-name=genomad     # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load fqtk/XXXX  ### Latest version is recommended.

   genomad end-to-end --cleanup --splits 8 genome.fa genomad_output /cluster/tufts/biocontainers/datasets/genomad_db/v1.7/
