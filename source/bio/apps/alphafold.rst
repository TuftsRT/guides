###########
 Alphafold
###########

**************
 Introduction
**************

``Alphafold`` is a protein structure prediction tool developed by
DeepMind (Google). It uses a novel machine learning approach to predict
3D protein structures from primary sequences alone. The source code is
available on Github_. It has been deployed in all RCAC clusters,
supporting both CPU and GPU.

It also relies on a huge database. The full database (~2.2TB) has been
downloaded and setup for users.

Protein struction prediction by alphafold is performed in the following
steps:

-  Search the amino acid sequence in uniref90 database by jackhmmer
   (using CPU)
-  Search the amino acid sequence in mgnify database by jackhmmer (using
   CPU)
-  Search the amino acid sequence in pdb70 database (for monomers) or
   pdb_seqres database (for multimers) by hhsearch (using CPU)
-  Search the amino acid sequence in bfd database and uniclust30
   (updated to uniref30 since v2.3.0) database by hhblits (using CPU)
-  Search structure templates in pdb_mmcif database (using CPU)
-  Search the amino acid sequence in uniprot database (for multimers) by
   jackhmmer (using CPU)
-  Predict 3D structure by machine learning (using CPU or GPU)
-  Structure optimisation with OpenMM (using CPU or GPU)

|  For more information, please check:
|  Home page: https://github.com/deepmind/alphafold

**********
 Versions
**********

-  2.3.0
-  2.3.1
-  2.3.2

**********
 Commands
**********

-  run_alphafold.sh

*******
 Usage
*******

The usage of Alphafold on our cluster is very straightford, users can
create a flagfile containing the database path information:

.. code::

   run_alphafold.sh --flagfile=full_db.ff --fasta_paths=XX --output_dir=XX ...

Users can check its detaied user guide in its Github_.

******************************************
 full_db_20230311.ff (for alphafold v2.3)
******************************************

Example contents of full_db_20231031.ff for monomer:

.. code::

   --db_preset=full_dbs
   --bfd_database_path=/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
   --data_dir=/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/
   --uniref90_database_path=/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/uniref90/uniref90.fasta
   --mgnify_database_path=/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/mgnify/mgy_clusters_2022_05.fa
   --uniref30_database_path=/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/uniref30/UniRef30_2021_03
   --pdb70_database_path=/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/pdb70/pdb70
   --template_mmcif_dir=/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/pdb_mmcif/mmcif_files
   --obsolete_pdbs_path=/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/pdb_mmcif/obsolete.dat
   --hhblits_binary_path=/usr/bin/hhblits
   --hhsearch_binary_path=/usr/bin/hhsearch
   --jackhmmer_binary_path=/usr/bin/jackhmmer
   --kalign_binary_path=/usr/bin/kalign

Example contents of full_db_20231031.ff for multimer:

.. code::

   --db_preset=full_dbs
   --bfd_database_path=/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
   --data_dir=/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/
   --uniref90_database_path=/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/uniref90/uniref90.fasta
   --mgnify_database_path=/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/mgnify/mgy_clusters_2022_05.fa
   --uniref30_database_path=/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/uniref30/UniRef30_2021_03
   --pdb_seqres_database_path=/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/pdb_seqres/pdb_seqres.txt
   --uniprot_database_path=/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/uniprot/uniprot.fasta
   --template_mmcif_dir=/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/pdb_mmcif/mmcif_files
   --obsolete_pdbs_path=/cluster/tufts/biocontainers/datasets/alphafold/db_20231031/pdb_mmcif/obsolete.dat
   --hhblits_binary_path=/usr/bin/hhblits
   --hhsearch_binary_path=/usr/bin/hhsearch
   --jackhmmer_binary_path=/usr/bin/jackhmmer
   --kalign_binary_path=/usr/bin/kalign

***********************
 Example job using CPU
***********************

To run alphafold using CPU:

.. code::

   #!/bin/bash
   #SBATCH -p PartitionName  # batch or your group's own partition
   #SBATCH -t 24:00:00
   #SBATCH -N 1
   #SBATCH -n 1
   #SBATCH -c 10
   #SBATCH --mem=64G
   #SBATCH --job-name=alphafold
   #SBATCH --mail-type=FAIL,BEGIN,END
   #SBATCH --error=%x-%J-%u.err
   #SBATCH --output=%x-%J-%u.out

   module purge
   module load alphafold/2.3.2

   run_alphafold.sh --flagfile=full_db_20231031.ff  \
       --fasta_paths=sample.fasta --max_template_date=2022-02-01 \
       --output_dir=af2_full_out --model_preset=monomer \
       --use_gpu_relax=False

***********************
 Example job using GPU
***********************

To run alphafold using GPU:

.. code::

   #!/bin/bash
   #SBATCH -p PartitionName  # gpu or preempt
   #SBATCH -t 24:00:00
   #SBATCH -N 1
   #SBATCH -n 1
   #SBATCH -c 10
   #SBATCH --mem=64G
   #SBATCH --gres=gpu:1
   #SBATCH --job-name=alphafold
   #SBATCH --mail-type=FAIL,BEGIN,END
   #SBATCH --error=%x-%J-%u.err
   #SBATCH --output=%x-%J-%u.out

   module purge
   module load alphafold/2.3.2

   run_alphafold.sh --flagfile=full_db_20231031.ff  \
       --fasta_paths=sample.fasta --max_template_date=2022-02-01 \
       --output_dir=af2_full_out --model_preset=monomer \
       --use_gpu_relax=True

.. _github: https://github.com/deepmind/alphafold/
