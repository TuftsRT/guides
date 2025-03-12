########
 Spades
########

**************
 Introduction
**************

Spades is an assembly toolkit containing various assembly pipelines.

**********
 Versions
**********

-  3.15.4
-  3.15.5

If you require a version newer than these two, please visit the `SPAdes
download and installation page
<https://ablab.github.io/spades/installation.html>`_

As of 2024-12-10, the latest version available is 4.0.0.

**********
 Commands
**********

-  coronaspades.py
-  metaplasmidspades.py
-  metaspades.py
-  metaviralspades.py
-  plasmidspades.py
-  rnaspades.py
-  rnaviralspades.py
-  spades-bwa
-  spades-convert-bin-to-fasta
-  spades-core
-  spades-corrector-core
-  spades-gbuilder
-  spades-gmapper
-  spades-gsimplifier
-  spades-hammer
-  spades_init.py
-  spades-ionhammer
-  spades-kmercount
-  spades-kmer-estimating
-  spades.py
-  spades-read-filter
-  spades-truseq-scfcorrection
-  truspades.py

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
   #SBATCH --job-name=spades      # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge
   module load spades/XXXX ### you can run *module avail spades* to check all available versions

********************************
 Using Multi-Threads for SPAdes
********************************

SPAdes is a genome/metagenome assembler designed to efficiently utilize
multiple CPU threads. Multi-threading significantly accelerates SPAdes
by parallelizing computationally intensive steps like k-mer
construction, error correction, and graph traversal.

How Multi-Threading Works in SPAdes:

-  **Parallelization**: SPAdes parallelizes tasks like k-mer generation,
   graph construction, and assembly across the threads specified with
   ``--threads``.

-  **Not Fully Parallel**: Some parts of SPAdes, such as graph
   simplification and certain I/O operations, are not fully
   parallelizable. As a result, doubling the number of threads may not
   always halve the runtime.

Best Practices for Using Multi-Threads with SPAdes:

-  **Determine Dataset Size**:
      -  Small datasets (e.g., bacterial genomes): 4–8 threads.
      -  Medium datasets (e.g., microbial communities): 8–16 threads.
      -  Large datasets (e.g., human or metagenomic assemblies): 16–32
         threads.

-  **Match SLURM Resources**:
      -  Ensure SLURM’s ``--cpus-per-task`` matches the ``--threads``
         argument in SPAdes.

-  **Monitor Memory Usage**:
      -  Check memory usage during assembly to ensure no resource
         wastage.

-  **Run Test Assemblies**:
      -  For new datasets, run SPAdes with fewer threads initially to
         estimate runtime and memory requirements.

.. code:: bash

   #!/bin/bash
   #SBATCH --job-name=SPAdes_Assembly
   #SBATCH --output=spades.%j.out
   #SBATCH --error=spades.%j.err
   #SBATCH --time=00-24:00:00  # Set appropriate runtime
   #SBATCH --cpus-per-task=16  # Match threads to SPAdes --threads
   #SBATCH --mem=64G           # Allocate sufficient memory
   #SBATCH -p batch            # Partition to use

   module load spades/3.15.5

   # Run SPAdes
   spades.py --meta \
          -1 /path/to/reads_R1.fastq.gz \
          -2 /path/to/reads_R2.fastq.gz \
          -o /path/to/output \
          --threads 16 \
          --memory 64

Using multi-threads is crucial for speeding up SPAdes, but it’s
important to balance thread count with available memory, I/O capacity,
and the complexity of your dataset to ensure efficient and successful
assemblies.

*****************
 Reference links
*****************

`SPAdes github repo <https://github.com/ablab/spades>`_

`SPAdes manual <https://ablab.github.io/spades/index.html>`_

`SPAdes paper <https://pmc.ncbi.nlm.nih.gov/articles/PMC3342519/>`_

`metaSPAdes paper: a new versatile metagenomic assembler
<https://pmc.ncbi.nlm.nih.gov/articles/PMC5411777/>`_
