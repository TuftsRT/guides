###########
 Metaphlan
###########

**************
 Introduction
**************

Metaphlan is computational tool for profiling the composition of
microbial communities (Bacteria, Archaea and Eukaryotes) from
metagenomic shotgun sequencing data (i.e. not 16S) with species-level.

**********
 Versions
**********

-  4.0.2
-  4.0.6

**********
 Commands
**********

-  add_metadata_tree.py
-  aggregate_scores_in_intervals.py
-  align_print_template.py
-  axt_extract_ranges.py
-  axt_to_fasta.py
-  axt_to_lav.py
-  axt_to_maf.py
-  bed_bigwig_profile.py
-  bed_build_windows.py
-  bed_complement.py
-  bed_count_by_interval.py
-  bed_count_overlapping.py
-  bed_coverage_by_interval.py
-  bed_coverage.py
-  bed_diff_basewise_summary.py
-  bed_extend_to.py
-  bed_intersect_basewise.py
-  bed_intersect.py
-  bed_merge_overlapping.py
-  bed_rand_intersect.py
-  bed_subtract_basewise.py
-  bnMapper.py
-  bowtie2
-  bowtie2-build
-  breadth_depth.py
-  cleanup-blastdb-volumes.py
-  consensus_aDNA.py
-  consensus.py
-  div_snp_table_chr.py
-  extract_markers.py
-  find_in_sorted_file.py
-  gene_fourfold_sites.py
-  get_scores_in_intervals.py
-  interval_count_intersections.py
-  interval_join.py
-  int_seqs_to_char_strings.py
-  lav_to_axt.py
-  lav_to_maf.py
-  line_select.py
-  lzop_build_offset_table.py
-  maf_build_index.py
-  maf_chop.py
-  maf_chunk.py
-  maf_col_counts_all.py
-  maf_col_counts.py
-  maf_count.py
-  maf_covered_ranges.py
-  maf_covered_regions.py
-  maf_div_sites.py
-  maf_drop_overlapping.py
-  maf_extract_chrom_ranges.py
-  maf_extract_ranges_indexed.py
-  maf_extract_ranges.py
-  maf_filter_max_wc.py
-  maf_filter.py
-  maf_gap_frequency.py
-  maf_gc_content.py
-  maf_interval_alignibility.py
-  maf_limit_to_species.py
-  maf_mapping_word_frequency.py
-  maf_mask_cpg.py
-  maf_mean_length_ungapped_piece.py
-  maf_percent_columns_matching.py
-  maf_percent_identity.py
-  maf_print_chroms.py
-  maf_print_scores.py
-  maf_randomize.py
-  maf_region_coverage_by_src.py
-  maf_select.py
-  maf_shuffle_columns.py
-  maf_species_in_all_files.py
-  maf_split_by_src.py
-  maf_thread_for_species.py
-  maf_tile_2bit.py
-  maf_tile_2.py
-  maf_tile.py
-  maf_to_axt.py
-  maf_to_concat_fasta.py
-  maf_to_fasta.py
-  maf_to_int_seqs.py
-  maf_translate_chars.py
-  maf_truncate.py
-  maf_word_frequency.py
-  mask_quality.py
-  merge_metaphlan_tables.py
-  metaphlan
-  mMK_bitset.py
-  nib_chrom_intervals_to_fasta.py
-  nib_intervals_to_fasta.py
-  nib_length.py
-  one_field_per_line.py
-  out_to_chain.py
-  plot_tree_graphlan.py
-  polymut.py
-  poly.py
-  prefix_lines.py
-  pretty_table.py
-  qv_to_bqv.py
-  random_lines.py
-  read_fastx.py
-  sample2markers.py
-  strainphlan
-  strain_transmission.py
-  sumlabels.py
-  sumtrees.py
-  table_add_column.py
-  table_filter.py
-  tfloc_summary.py
-  ucsc_gene_table_to_intervals.py
-  wiggle_to_array_tree.py
-  wiggle_to_binned_array.py
-  wiggle_to_chr_binned_array.py
-  wiggle_to_simple.py

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
   #SBATCH --job-name=metaphlan   # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load metaphlan/XXXX     ### Latest version is recommended.
