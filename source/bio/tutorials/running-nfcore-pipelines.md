---
tags: bioinformatics
---

# Running nf-core pipelines on Tufts HPC

Shirley Li, Bioinformatician, TTS Research Technology
xue.li37@tufts.edu

Date: 2026-02-25

You have two options for running the nf-core pipeline on Tufts HPC: through Open OnDemand or the command-line interface.

## 1. Open Ondemand

Navigate to [Open Ondemand](https://ondemand-p01.pax.tufts.edu/), under `nf-core pipelines`, choose the appropriate application for the pipeline you wish to execute.

## 2. Command Line Interface

Load the necessary modules

```
module load nextflow
module load singularity/4.3.4
```

Run the pipeline with this command:

```
nextflow run nf-core/rnaseq -profile tufts ...
```

For a comprehensive list of available pipelines, visit the [nf-core website](https://nf-co.re/pipelines). There are currently more than 100 pipelines available through nf-core.

### Example slurm script

```
#!/bin/bash
#SBATCH --job-name=rnaseq_nfcore
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=rnaseq_%j.log
#SBATCH --error=rnaseq_%j.err

echo "Starting job at $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"

module load nextflow
module load singularity/4.3.4

genome=/cluster/tufts/workshop/public/2026spring/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gtf=/cluster/tufts/workshop/public/2026spring/reference/Homo_sapiens.GRCh38.111.gtf.gz

nextflow run nf-core/rnaseq \
    -r 3.22.2 \
    -profile tufts \
    --input samplesheet.csv \
    --fasta $genome \
    --gtf $gtf \
    --star_index /cluster/tufts/workshop/public/2026spring/star_index/  \
    --salmon_index /cluster/tufts/workshop/public/2026spring/salmon_index/
    --outdir rnaseq


echo "Finished at $(date)"
```

## 3. Key Results

If you define your output directory as:

```
--outdir rnaseq
```

then all pipeline results will be written inside:

```
rnaseq/
```

This directory will contain subfolders such as:

```
rnaseq/
├── multiqc/
├── star_salmon/
├── pipeline_info/
├── fastqc/
├── trimgalore/
└── work/   (intermediate files, DELETE after COMPLETION to save space)
```

### 3.1 Summary Report (Most Important)

```
multiqc/
multiqc/star_salmon/multiqc_report.html
```

Open this first.

### 3.2 Gene-Level Counts

```
salmon.merged.gene_counts.tsv   → gene-level counts (used for differential expression)
```

```
salmon.merged.gene_tpm.tsv      → normalized expression (for visualization)
```

### 3.3 Pipeline Metadata (Reproducibility)

```
pipeline_info/nf_core_rnaseq_software_mqc_versions.yml
```

Contains software versions

#### Open OnDemand App

In your working directory, this file list all parameter used

```
nf-params.json
```

## 4. Logs and Troubleshooting

### 4.1 Main Nextflow Log

```
.nextflow.log
```

Monitor during run:

```
tail -f .nextflow.log
```

It locates in the workding dir, not the outdir.

### 4.2 SLURM Log Files

```
rnaseq_<jobid>.log   # check progress
rnaseq_<jobid>.err   # check errors
```

#### Open OnDemand App

In session ID folder:

```
output.log  # check progress
```

### 4.3 Per-Process Logs (Advanced)

Located inside:

```
work/<hash>/
```

Files:

```
.command.log
.command.err
.command.out
```
