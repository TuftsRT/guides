---
tags: bioinformatics
---

# Getting started with bioinformatics tools

Shirley Li, Bioinformatician, TTS Research Technology
xue.li37@tufts.edu

Date: 2024-11-01

Warning: This document serves as a guide and is not intended to be the definitive standard for bioinformatics analysis.

For your specific research needs, always consult the original references and adhere to established best practices. Additionally, ensure you visit the official websites of the tools listed herein to confirm their usage instructions and obtain the most accurate and up-to-date information.

Feel free to reach out to tts-research@tufts.edu, if you notice any errors in this documentation.

## Introduction

This guide provides an overview of the bioinformatics tools currently available on our cluster, organized by specific research areas. Whether you're analyzing RNA-Seq data, exploring whole genome sequencing, conducting ChIP-seq experiments, delving into meta-genomics, or predicting protein structures, this resource aims to help you identify the most suitable tools for your research needs. Each section lists tools along with their versions to facilitate easy module loading and integration into your workflows.

**Note: Many tools are integrated into interactive apps available on [Open Ondemand](https://ondemand.pax.tufts.edu)**\
**Note: For standardized and reproducible workflows, we recommend using nf-core pipelines. Instructions and guidelines are available here (link to nf-core tutorial will be inserted)**

## RNA-Seq Data Analysis

- **QC**: fastqc/0.11.8, trim-galore/0.6.10, multiqc/1.7.0
- **Alignment**: STAR/2.7.10a, bowtie2/2.5.1, hisat2/2.2.1
- **Quantification**: cufflinks/2.2.1, rsem/1.3.3, salmon/1.10.1
- **Differential Expression Analysis**: DESeq2 (used via R packages), edgeR (used via R packages)
- **Visualization**: R/4.3.0

## ChIP-seq Analysis

- **QC**: fastqc/0.11.8, trim-galore/0.6.10
- **Alignment**: bowtie2/2.5.1, STAR/2.7.10a
- **Peak Calling**: macs2/2.2.7.1
- **Motif Analysis**: homer/4.11
- **Visualization**: R/4.3.0

## Whole Genome Sequencing (WGS) Analysis

- **QC**: fastqc/0.11.8, trim-galore/0.6.10, multiqc/1.7.0
- **Alignment**: bwa/0.7.17, bowtie2/2.5.1
- **Variant Calling and Management**: bcftools/1.14, freebayes/1.1.0, Picard/2.26.10, GATK/4.2.6.1, vcftools/0.1.16
- **Structural Variant Analysis**: delly/0.8.1 (not listed, commonly used)
- **Annotation**: snpeff/4.3t (not listed, commonly used)
- **Visualization**: IGV (Open Ondemand or Local)
- **Accelerated Genomic Processing with GPU**: parabricks/4.2.1-1

## Genetic Association Studies

- **Genotype and Phenotype Analysis**: plink/1.90b6.21, gcta/1.25.2

## Metagenomics

- **QC**: fastqc/0.11.8, multiqc/1.7.0
- **Assembly**: megahit/1.2.9
- **Taxonomic Classification**: kraken2/2.1.1 (not listed, commonly used)
- **Microbial Community Analysis**: mothur/1.48.0, QIIME/2023.9
- **Functional Annotation**: humann/3.8
- **Protein Sequence Alignment**: diamond/2.1.6
- **Visualization**: R/4.3.0
- **Workflow Management**: biobakery_workflows/3.0.0.a.7

## Protein Structure Prediction

- **Structure Modeling**: alphafold/2.3.2, rosetta/3.13,
- **Homology Modeling**: modeller/9.25 (not listed, commonly used)
- **Molecular Dynamics Simulations**: gromacs/2023.3_gcc_9.3.0
- **Visualization**: pymol (Open Ondemand or Local)

## Alignment Tools

- **DNA Sequence Alignment**: bwa/0.7.17, bowtie2/2.5.1
- **RNA-Seq Alignment**: hisat2/2.2.1, STAR/2.7.10a
- **Short Read Alignment**: bowtie2/2.5.1
- **Long Read Alignment**: minimap2/2.15

## Ancient DNA Analysis

- **Read Alignment**: bwa/0.7.17, bowtie2/2.5.1
- **Data Authentication and Decontamination**: angsd/0.940, mapDamage/2.0.9 (not listed)
- **Genotype Calling and Analysis**: samtools/1.9, GATK/4.2.6.1
- **Population Genetics**: plink/1.90b6.21, ADMIXTURE/1.3.0 (not listed)

## Nextflow and nf-core Tools

- **Workflow management**: nextflow/23.10.0
- **Standardized Pipelines**: nf-core/2.10
- **Containerization Tool**: singularity/3.8.4 (default)
- **nf-core Modules**: nf-core-ampliseq/2.8.0, nf-core-ampliseq/2.9.0, nf-core-atacseq/2.1.2, nf-core-bacass/2.2.0, nf-core-bamtofastq/2.1.1, nf-core-chipseq/2.0.0, nf-core-differentialabundance/1.4.0, nf-core-differentialabundance/1.5.0, nf-core-eager/2.5.1, nf-core-fetchngs/1.11.0, nf-core-fetchngs/1.12.0, nf-core-funcscan/1.1.4, nf-core-hic/2.1.0, nf-core-mag/2.5.2, nf-core-mag/2.5.4, nf-core-mag/3.0.0, nf-core-metatdenovo/1.0.0, nf-core-methylseq/2.6.0, nf-core-nanoseq/3.1.0, nf-core-nanostring/1.2.1, nf-core-pangenome/1.1.0, nf-core-pangenome/1.1.1, nf-core-raredisease/2.0.1, nf-core-rnafusion/3.0.1, nf-core-rnaseq/3.14.0, nf-core-rnasplice/1.0.2, nf-core-rnasplice/1.0.3, nf-core-rnasplice/1.0.4, nf-core-sarek/3.4.0, nf-core-sarek/3.4.1, nf-core-scrnaseq/2.5.1, nf-core-smrnaseq/2.3.0, nf-core-taxprofiler/1.1.5, nf-core-taxprofiler/1.1.6, nf-core-taxprofiler/1.1.7, nf-core-viralrecon/2.6.0.
