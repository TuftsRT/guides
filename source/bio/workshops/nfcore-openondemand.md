# Running nf-core RNA-seq on Open OnDemand and CLI

2026-02-20

Shirley Li: xue.li37@tufts.edu



Pipeline version used in workshop: [**nf-core/rnaseq 3.22.2**](https://nf-co.re/rnaseq/3.22.2/)

# 1. Prepare Input Files

## 1.1 FASTQ Files

Place all FASTQ files in your working directory.

Example (paired-end):

```
SRX1693951_SRR3362661_1_sub.fastq.gz
SRX1693951_SRR3362661_2_sub.fastq.gz
...
SRX1693956_SRR3362666_2_sub.fastq.gz
```

Our sample fastq file located here:

```
/cluster/tufts/workshop/public/2026spring/nfcore/fastq/
```





## 1.2 Sample Sheet (`samplesheet.csv`)

Required format:

```
sample,fastq_1,fastq_2,strandedness
GFPkd_1,SRX1693951_SRR3362661_1_sub.fastq.gz,SRX1693951_SRR3362661_2_sub.fastq.gz,auto
GFPkd_2,SRX1693952_SRR3362662_1_sub.fastq.gz,SRX1693952_SRR3362662_2_sub.fastq.gz,auto
GFPkd_3,SRX1693953_SRR3362663_1_sub.fastq.gz,SRX1693953_SRR3362663_2_sub.fastq.gz,auto
PRMT5kd_1,SRX1693954_SRR3362664_1_sub.fastq.gz,SRX1693954_SRR3362664_2_sub.fastq.gz,auto
PRMT5kd_2,SRX1693955_SRR3362665_1_sub.fastq.gz,SRX1693955_SRR3362665_2_sub.fastq.gz,auto
PRMT5kd_3,SRX1693956_SRR3362666_1_sub.fastq.gz,SRX1693956_SRR3362666_2_sub.fastq.gz,auto
```

Notes:

- Paired-end files must match correctly.
- Paths can be relative or absolute (recommended).



On HPC:

```
/cluster/tufts/workshop/public/2026spring/nfcore/samplesheet.csv
```



## 1.3 Reference Files

You must provide:

- Genome FASTA file
- Gene annotation GTF file

Reference files can be either **remote (URL)** or **local (file path on the cluster)**.

### Remote URLs

Example (Ensembl GRCh38 release 111):

FASTA (Genome):

```
https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

GTF (Annotation):

```
https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz
```

### Local example (For workshop use only). 

**This directory may be deleted after the workshop. Keep your own copy if you plan to rerun the analysis.**

FASTA (Genome):

```
/cluster/tufts/workshop/public/2026spring/star_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

GTF (Annotation):

```
/cluster/tufts/workshop/public/2026spring/star_index/Homo_sapiens.GRCh38.111.gtf
```





# 2. Running via Open OnDemand

Log in to:

Open OnDemand → nf-core pipelines → rnaseq (version 3.22.2)

Configure:

- Input: `samplesheet.csv`
- FASTA: (URL or local file path listed above)
- GTF: (URL or local file path listed above)
- Outdir: `rnaseq`
- Working directory: your project directory
- kraken_db: `/cluster/tufts/biocontainers/datasets/kraken2/k2_standard_20251015`
- star_index: `/cluster/tufts/workshop/public/2026spring/star_index/ `
- salmon_index: `/cluster/tufts/workshop/public/2026spring/salmon_index/ `
- skip_pseudo_alignment: true

Other parameters can remain default unless discussed in workshop.

Submit job.



# 3. Running via CLI (SLURM)

Create a SLURM script (example: `run_rnaseq.slurm`)

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

# -resume 

# --star_index is provided for the workshop to avoid rebuilding the index and save time.
# You do not need this flag when running the pipeline independently.

echo "Finished at $(date)"
```

Submit job:

```
sbatch run_rnaseq.slurm
```



# 4. Key Results

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

## 4.1 Summary Report (Most Important)

```
multiqc/
multiqc/star_salmon/multiqc_report.html 
```

Open this first.

## 4.2 Gene-Level Counts

```
salmon.merged.gene_counts.tsv   → gene-level counts (used for differential expression)
```

```
salmon.merged.gene_tpm.tsv      → normalized expression (for visualization)
```



## 4.3 Pipeline Metadata (Reproducibility)

```
pipeline_info/nf_core_rnaseq_software_mqc_versions.yml
```

Contains software versions 



### Open OnDemand App

In your working directory, this file list all parameter used

```
nf-params.json
```





# 5. Logs and Troubleshooting

## 5.1 Main Nextflow Log

```
.nextflow.log
```

Monitor during run:

```
tail -f .nextflow.log
```

It locates in the workding dir, not the outdir. 

## 5.2 SLURM Log Files

```
rnaseq_<jobid>.log   # check progress
rnaseq_<jobid>.err   # check errors
```

### Open OnDemand App

In session ID folder:

```
output.log  # check progress
```



## 5.3 Per-Process Logs (Advanced)

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



# 6. Running With Your Own Data

To run your own dataset, you need:

- FASTQ files
- A properly formatted `samplesheet.csv`
- FASTA reference genome
- GTF annotation file

Steps:

1. Upload FASTQ files to your working directory.
2. Create a `samplesheet.csv`.
3. Adjust reference paths if needed.
4. Launch the pipeline via Open OnDemand or CLI.