####################
 Pepper_deepvariant
####################

**************
 Introduction
**************

PEPPER is a genome inference module based on recurrent neural networks
that enables long-read variant calling and nanopore assembly polishing
in the PEPPER-Margin-DeepVariant pipeline.

**********
 Versions
**********

-  r0.8
-  r0.8-gpu

**********
 Commands
**********

-  call_variants
-  freeze_graph
-  make_examples
-  model_eval
-  model_train
-  multisample_make_examples
-  pepper
-  pepper_train
-  pepper_variant
-  pepper_variant_train
-  postprocess_variants
-  run_deepvariant
-  run_deepvariant.py
-  run_pepper_margin_deepvariant
-  run-prereq.sh
-  runtime_by_region_vis
-  settings.sh
-  show_examples
-  vcf_stats_report

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
   #SBATCH --job-name=pepper_deepvariant  # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load pepper_deepvariant/XXXX    ### Latest version is recommended.
