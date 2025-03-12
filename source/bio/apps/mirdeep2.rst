##########
 Mirdeep2
##########

**************
 Introduction
**************

miRDeep2 discovers active known or novel miRNAs from deep sequencing
data (Solexa/Illumina, 454, ...).

**********
 Versions
**********

-  2.0.1.3

**********
 Commands
**********

-  bwa_sam_converter.pl
-  clip_adapters.pl
-  collapse_reads_md.pl
-  convert_bowtie_output.pl
-  excise_precursors_iterative_final.pl
-  excise_precursors.pl
-  extract_miRNAs.pl
-  fastaparse.pl
-  fastaselect.pl
-  fastq2fasta.pl
-  find_read_count.pl
-  geo2fasta.pl
-  get_mirdeep2_precursors.pl
-  illumina_to_fasta.pl
-  make_html2.pl
-  make_html.pl
-  mapper.pl
-  mirdeep2bed.pl
-  miRDeep2_core_algorithm.pl
-  miRDeep2.pl
-  parse_mappings.pl
-  perform_controls.pl
-  permute_structure.pl
-  prepare_signature.pl
-  quantifier.pl
-  remove_white_space_in_id.pl
-  rna2dna.pl
-  samFLAGinfo.pl
-  sam_reads_collapse.pl
-  sanity_check_genome.pl
-  sanity_check_mapping_file.pl
-  sanity_check_mature_ref.pl
-  sanity_check_reads_ready_file.pl
-  select_for_randfold.pl
-  survey.pl

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
   #SBATCH --job-name=mirdeep2    # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load mirdeep2/XXXX      ### Latest version is recommended.
