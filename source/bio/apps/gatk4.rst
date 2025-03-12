#######
 Gatk4
#######

**************
 Introduction
**************

`GATK (Genome Analysis Toolkit)
<https://gatk.broadinstitute.org/hc/en-us>`_ is a collection of
command-line tools for analyzing high-throughput sequencing data with a
primary focus on variant discovery.

**********
 Versions
**********

-  4.2.6.1
-  4.3.0.0
-  4.5.0.0

**********
 Commands
**********

-  gatk

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
   #SBATCH --job-name=gatk4       # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load gatk4/XXXX ### Latest version is recommended.


   # Define input and output files
   REF_GENOME="reference.fasta"
   INPUT_BAM="sample.bam"
   OUTPUT_VCF="sample.vcf"

   # Run GATK HaplotypeCaller
   gatk HaplotypeCaller \
     -R $REF_GENOME \
     -I $INPUT_BAM \
     -O $OUTPUT_VCF \
     --native-pair-hmm-threads 4

Best Practices:

-  Always test your scripts with a small dataset before scaling up.
-  Optimize memory and CPU usage based on your data size and cluster
   configuration.
-  Monitor job performance using `squeue` or similar cluster tools.

************
 References
************

`Tool Documentation Index for v4.5.0.0
<https://gatk.broadinstitute.org/hc/en-us/articles/21904996835867--Tool-Documentation-Index>`_

`Human genome reference builds - GRCh38 or hg38 - b37 - hg19
<https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19>`_

`GATK Community Forum
<https://gatk.broadinstitute.org/hc/en-us/community/topics>`_

`GATK citations
<https://gatk.broadinstitute.org/hc/en-us/articles/360035530852-How-should-I-cite-GATK-in-my-own-publications>`_
