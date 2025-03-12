########
 Bedops
########

**************
 Introduction
**************

Bedops is a software package for manipulating and analyzing genomic
interval data.

**********
 Versions
**********

-  2.4.39

**********
 Commands
**********

-  bam2bed
-  bam2bed-float128
-  bam2bed_gnuParallel
-  bam2bed_gnuParallel-float128
-  bam2bed_gnuParallel-megarow
-  bam2bed_gnuParallel-typical
-  bam2bed-megarow
-  bam2bed_sge
-  bam2bed_sge-float128
-  bam2bed_sge-megarow
-  bam2bed_sge-typical
-  bam2bed_slurm
-  bam2bed_slurm-float128
-  bam2bed_slurm-megarow
-  bam2bed_slurm-typical
-  bam2bed-typical
-  bam2starch
-  bam2starch-float128
-  bam2starch_gnuParallel
-  bam2starch_gnuParallel-float128
-  bam2starch_gnuParallel-megarow
-  bam2starch_gnuParallel-typical
-  bam2starch-megarow
-  bam2starch_sge
-  bam2starch_sge-float128
-  bam2starch_sge-megarow
-  bam2starch_sge-typical
-  bam2starch_slurm
-  bam2starch_slurm-float128
-  bam2starch_slurm-megarow
-  bam2starch_slurm-typical
-  bam2starch-typical
-  bedextract
-  bedextract-float128
-  bedextract-megarow
-  bedextract-typical
-  bedmap
-  bedmap-float128
-  bedmap-megarow
-  bedmap-typical
-  bedops
-  bedops-float128
-  bedops-megarow
-  bedops-typical
-  closest-features
-  closest-features-float128
-  closest-features-megarow
-  closest-features-typical
-  convert2bed
-  convert2bed-float128
-  convert2bed-megarow
-  convert2bed-typical
-  gff2bed
-  gff2bed-float128
-  gff2bed-megarow
-  gff2bed-typical
-  gff2starch
-  gff2starch-float128
-  gff2starch-megarow
-  gff2starch-typical
-  gtf2bed
-  gtf2bed-float128
-  gtf2bed-megarow
-  gtf2bed-typical
-  gtf2starch
-  gtf2starch-float128
-  gtf2starch-megarow
-  gtf2starch-typical
-  gvf2bed
-  gvf2bed-float128
-  gvf2bed-megarow
-  gvf2bed-typical
-  gvf2starch
-  gvf2starch-float128
-  gvf2starch-megarow
-  gvf2starch-typical
-  psl2bed
-  psl2bed-float128
-  psl2bed-megarow
-  psl2bed-typical
-  psl2starch
-  psl2starch-float128
-  psl2starch-megarow
-  psl2starch-typical
-  rmsk2bed
-  rmsk2bed-float128
-  rmsk2bed-megarow
-  rmsk2bed-typical
-  rmsk2starch
-  rmsk2starch-float128
-  rmsk2starch-megarow
-  rmsk2starch-typical
-  sam2bed
-  sam2bed-float128
-  sam2bed-megarow
-  sam2bed-typical
-  sam2starch
-  sam2starch-float128
-  sam2starch-megarow
-  sam2starch-typical
-  sort-bed
-  sort-bed-float128
-  sort-bed-megarow
-  sort-bed-typical
-  starch
-  starchcat
-  starchcat-float128
-  starchcat-megarow
-  starchcat-typical
-  starchcluster_gnuParallel
-  starchcluster_gnuParallel-float128
-  starchcluster_gnuParallel-megarow
-  starchcluster_gnuParallel-typical
-  starchcluster_sge
-  starchcluster_sge-float128
-  starchcluster_sge-megarow
-  starchcluster_sge-typical
-  starchcluster_slurm
-  starchcluster_slurm-float128
-  starchcluster_slurm-megarow
-  starchcluster_slurm-typical
-  starch-diff
-  starch-diff-float128
-  starch-diff-megarow
-  starch-diff-typical
-  starch-float128
-  starch-megarow
-  starchstrip
-  starchstrip-float128
-  starchstrip-megarow
-  starchstrip-typical
-  starch-typical
-  switch-BEDOPS-binary-type
-  unstarch
-  unstarch-float128
-  unstarch-megarow
-  unstarch-typical
-  update-sort-bed-migrate-candidates
-  update-sort-bed-migrate-candidates-float128
-  update-sort-bed-migrate-candidates-megarow
-  update-sort-bed-migrate-candidates-typical
-  update-sort-bed-slurm
-  update-sort-bed-slurm-float128
-  update-sort-bed-slurm-megarow
-  update-sort-bed-slurm-typical
-  update-sort-bed-starch-slurm
-  update-sort-bed-starch-slurm-float128
-  update-sort-bed-starch-slurm-megarow
-  update-sort-bed-starch-slurm-typical
-  vcf2bed
-  vcf2bed-float128
-  vcf2bed-megarow
-  vcf2bed-typical
-  vcf2starch
-  vcf2starch-float128
-  vcf2starch-megarow
-  vcf2starch-typical
-  wig2bed
-  wig2bed-float128
-  wig2bed-megarow
-  wig2bed-typical
-  wig2starch
-  wig2starch-float128
-  wig2starch-megarow
-  wig2starch-typical

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
   #SBATCH --job-name=bedops      # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load bedops/XXXX        ### Latest version is recommended.
