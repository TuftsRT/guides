########
 Humann
########

**************
 Introduction
**************

Humann is a pipeline for efficiently and accurately profiling the
presence/absence and abundance of microbial pathways in a community from
metagenomic or metatranscriptomic sequencing data (typically millions of
short DNA/RNA reads).

**********
 Versions
**********

-  3.8

**********
 Commands
**********

-  humann
-  humann3
-  humann3_databases
-  humann_associate
-  humann_barplot
-  humann_benchmark
-  humann_build_custom_database
-  humann_config
-  humann_databases
-  humann_genefamilies_genus_level
-  humann_humann1_kegg
-  humann_infer_taxonomy
-  humann_join_tables
-  humann_reduce_table
-  humann_regroup_table
-  humann_rename_table
-  humann_renorm_table
-  humann_rna_dna_norm
-  humann_split_stratified_table
-  humann_split_table
-  humann_strain_profiler
-  humann_test
-  humann_unpack_pathways

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
   #SBATCH --job-name=humann      # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load humann/XXXX        ### Latest version is recommended.
