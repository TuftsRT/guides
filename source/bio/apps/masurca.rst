#########
 Masurca
#########

**************
 Introduction
**************

The MaSuRCA (Maryland Super Read Cabog Assembler) genome assembly and
analysis toolkit contains of MaSuRCA genome assembler, QuORUM error
corrector for Illumina data, POLCA genome polishing software, Chromosome
scaffolder, jellyfish mer counter, and MUMmer aligner.

**********
 Versions
**********

-  4.0.9
-  4.1.0

**********
 Commands
**********

-  build_human_reference.sh
-  chromosome_scaffolder.sh
-  close_gaps.sh
-  close_scaffold_gaps.sh
-  correct_with_k_unitigs.sh
-  deduplicate_contigs.sh
-  deduplicate_unitigs.sh
-  eugene.sh
-  extract_chrM.sh
-  filter_library.sh
-  final_polish.sh
-  fix_unitigs.sh
-  fragScaff.sh
-  masurca
-  mega_reads_assemble_cluster2.sh
-  mega_reads_assemble_cluster.sh
-  mega_reads_assemble_polish.sh
-  mega_reads_assemble_ref.sh
-  parallel_delta-filter.sh
-  polca.sh
-  polish_with_illumina_assembly.sh
-  recompute_astat_superreads_CA8.sh
-  recompute_astat_superreads.sh
-  reconcile_alignments.sh
-  refine.sh
-  resolve_trio.sh
-  run_ECR.sh
-  samba.sh
-  splitScaffoldsAtNs.sh

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
   #SBATCH --job-name=masurca     # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load masurca/XXXX       ### Latest version is recommended.
