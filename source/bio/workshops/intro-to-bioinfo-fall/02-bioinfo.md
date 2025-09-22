# Introduction to Bioinformatics on Tufts HPC

Author: Shirley Li, xue.li37@tufts.edu

Date: 2025-02

## Prerequisites

- Basic understanding of biology and bioinformatics

  - Familiarity with common bioinformatics data formats, including [FASTQ](https://knowledge.illumina.com/software/general/software-general-reference_material-list/000002211), [FASTA](https://pacbiofileformats.readthedocs.io/en/11.0/FASTA.html), [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf), [BAM](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm), and [GTF files](https://useast.ensembl.org/info/website/upload/gff.html).

- [Familiarity with the command line](https://tuftsdatalab.github.io/tuftsWorkshops/2024_workshops/2024_bioinformatics201/linux/00_overview/)

- [Access to an HPC cluster](https://it.tufts.edu/researchtechnology.tufts.edu) (e.g., login credentials, necessary software installations)

## Bioinformatics modules

### On the cluster

Use `module avail` to check the full list of tools available on the cluster. Below are some commonly used tools:

```
   abcreg/0.1.0                         kallisto/0.48.0                     (D)    orthofinder/2.5.5
   abyss/2.3.7                          kneaddata/0.12.0                           pandaseq/2.11
   alphafold/2.3.0                      kraken2/2.1.3                              parabricks/4.0.0-1
   alphafold/2.3.1                      krakentools/1.2                            parabricks/4.2.1-1
   alphafold/2.3.2                      macs2/2.2.7.1
   amplify/2.0.0                        macs3/3.0.0a6                              pepper_deepvariant/r0.8
   angsd/0.939                          masurca/4.0.9                              petitefinder/cpu
   angsd/0.940                   (D)    masurca/4.1.0                       (D)    picard/2.25.1
   bakta/1.9.3                          medaka/1.11.1                              picard/2.26.10
   bbmap/38.93                          megahit/1.2.9                              plink/1.90b6.21
   bbmap/38.96                   (D)    meme/5.5.5                                 plink2/2.00a2.3
   bbtools/39.00                        metaphlan/4.0.2                            polypolish/0.5.0
   bcftools/1.13                        metaphlan/4.0.6                     (D)    preseq/3.2.0
   bcftools/1.14                        miniasm/0.3_r179                           prokka/1.14.6
   bcftools/1.17                        minimap2/2.26                       (D)    qiime2/2023.2
   bcftools/1.20                 (D)    minipolish/0.1.3                           qiime2/2023.5
   beast2/2.6.3                         mirdeep2/2.0.1.3                           qiime2/2023.7
   beast2/2.6.4                         mirge3/0.1.4                               qiime2/2023.9
   beast2/2.6.6                  (D)    mothur/1.46.0                              qiime2/2024.2
   ... ...

```

### A few tips

1. Before installing your own tools, check if they are already available on the cluster using the `module avail` command.

1. Always be aware of the software versions, especially when using scripts from colleagues.

1. For less common tools, consider installing them yourself to ensure you have full control over the version and availability.

   **If you need to install a less commonly used tool, it's best to handle the installation yourself to ensure proper maintenance. Follow [this tutorial](https://tuftsdatalab.github.io/tuftsWorkshops/2024_workshops/2024_bioinformatics301/01_source/) to install your own tool**

## Using the Open OnDemand App

You can access Open OnDemand through [this link](https://ondemand.pax.tufts.edu/). 

If you have access to the new cluster, you can access the [new Open OnDemand server](http://ondemand-prod.pax.tufts.edu)


### RStudio and Apps

RStudio Pax, use R/4.4.2 which has the most comprehensive packages installed (1300+).

**How to initiate an R job**

1. Log in to Open Ondemand.
1. Navigate to `interactive apps` and select `RStudio Pax`
1. Specify the required resources:
   - **Number of hours**
   - **Number of CPU cores**
   - **Amount of memory**
   - **CPU partition** (set to `batch`)
   - **R version** (latest available: 4.4.2)
1. Click `Launch` to submit your job to the queue.
1. Wait a few minutes until your job starts running.
1. Click `Connect to Rstudio server`
1. In RStudio, go to the `Packages` tab on the right to check the installed packages.

**Installing R Packages**

Refer to our [previous workshop materials](https://tuftsdatalab.github.io/tuftsWorkshops/2024_workshops/2024_bioinformatics301/02_Rpackage/) for detailed instructions on installing R packages.

### Other Apps

We also provide other applications like `Jupyter Bioinfo`, `JupyterLab`, `Jupyter Notebook`, `IGV`, and `Galaxy` to support your daily research activities.

### [nf-core pipelines](https://nf-co.re/pipelines)

- On cluster, use `module avail nf-core` to get the list of nf-core pipelines deployed on cluster.

- On Open OnDemand, you can go to `bioinformatics apps` to find out what has been installed.

- [Quick Start Guide to Using the nf-core Pipeline](https://tuftsrt.github.io/guides/dev/bio/tutorials/doc/Running_nfcore_pipelines_at_Tufts_HPC.html)

## Writing Bioinformatics job script

Let’s explore how to write a SLURM script, using STAR alignment as an example.

### 1. Prepare the SLURM Script

Create a file named `run_star.sh` and add the following content:

```
#!/bin/bash
#SBATCH -J STAR_JOB             # Job name
#SBATCH --time=12:00:00         # Maximum runtime (D-HH:MM:SS format)
#SBATCH -p batch                # Partition (queue) to submit the job to
#SBATCH -n 1                    # Number of tasks (1 task in this case)
#SBATCH --mem=32g               # Memory allocation (32 GB)
#SBATCH --cpus-per-task=8       # Number of CPU cores allocated for the task
#SBATCH --output=STAR.%j.out    # Standard output file (%j = Job ID)
#SBATCH --error=STAR.%j.err     # Standard error file (%j = Job ID)
#SBATCH --mail-type=ALL         # Notifications for job status (start, end, fail)
#SBATCH --mail-user=utln@tufts.edu  # Your email address for notifications

# Load necessary module
module load star/2.7.11b

# Create output directory
mkdir -p star_output

# Run STAR alignment for single-end reads
STAR --genomeDir ./reference_data/reference_index/ \
     --readFilesIn ./raw_fastq/Irrel_kd_1.subset.fq \
     --outFileNamePrefix ./star_output/ \
     --runThreadN 8
```

### 2. Submit the Job Script

If you're running the script **directly in the terminal**, you need to make it executable first:

```
chmod +x run_star.sh
```

However, **SLURM does not require execution permissions**, so you can submit the job as-is using:

```
sbatch run_star.sh
```

### 3. Monitor Job Status

Use the following command to check the job status:

```
squeue -u yourusername
```

### STAR Workflow Details

#### Loading Modules

Before using STAR, load the appropriate module:

```
module load star/2.7.11b
```

#### Generating the Genome Index

Before aligning reads, generate the genome index using a reference genome (`.fa`) and an annotation file (`.gtf`):

```
STAR --runMode genomeGenerate \
     --genomeDir ./reference_data/ \
     --genomeFastaFiles ./reference_data/chr1.fa \
     --sjdbGTFfile ./reference_data/chr1-hg19_genes.gtf \
     --runThreadN 8
```

#### Aligning Reads

##### For Single-End Reads:

```
STAR --genomeDir ./reference_data/ \
     --readFilesIn ./raw_fastq/Irrel_kd_1.subset.fq \
     --outFileNamePrefix ./star_output/ \
     --runThreadN 8
```

##### For Paired-End Reads:

```
STAR --genomeDir ./reference_data/ \
     --readFilesIn /path/to/read1.fastq /path/to/read2.fastq \
     --outFileNamePrefix ./star_output/ \
     --runThreadN 8
```

#### Output Files

After running STAR, the output directory will contain:

```
Aligned.out.sam  Log.final.out  Log.out  Log.progress.out  SJ.out.tab
```

- **`Aligned.out.sam`**: Contains alignment data.
- **`Log.final.out`**: Summarizes alignment metrics.

#### Additional Tips

- Always test commands interactively before incorporating them into job scripts.

- Use the SLURM `--time`, `--mem`, and `--cpus-per-task` options to optimize resource allocation.

- Check the SLURM output and error files for troubleshooting.

### Run job with GPU node

#### Interactive session

```
srun -p preempt -n 1 --time=04:00:00 --mem=20G --gres=gpu:1 --pty /bin/bash
```

You can also specify which gpu node you would like to run jobs on

```
srun -p preempt -n 1 --time=04:00:00 --mem=20G --gres=gpu:a100:1 --pty /bin/bash
```

### Submit jobs to queue

Example script: `align.sh` using parabricks to do the alignment.

```
#!/bin/bash
#SBATCH -J fq2bam_alignment          # Job name
#SBATCH -p preempt                   # Submit to the 'preempt' partition (modify based on your cluster setup)
#SBATCH --gres=gpu:1                 # Request 1 GPU for accelerated processing
#SBATCH -n 2                         # Number of tasks (2 in this case)
#SBATCH --mem=60g                    # Memory allocation (60GB)
#SBATCH --time=02:00:00              # Maximum job run time (2 hours)
#SBATCH --cpus-per-task=20           # Number of CPU cores allocated per task
#SBATCH --output=alignment.%j.out    # Standard output file (with job ID %j)
#SBATCH --error=alignment.%j.err     # Standard error file (with job ID %j)
#SBATCH --mail-type=ALL              # Email notifications for all job states (begin, end, fail)
#SBATCH --mail-user=utln@tufts.edu   # Email address for notifications

# Load necessary modules
nvidia-smi                              # Show GPU information (optional for logging)
module load parabricks/4.0.0-1          # Load Parabricks module for GPU-accelerated alignment

# Define variables
genome_reference="/path/to/reference_genome"      # Path to the reference genome (.fasta)
input_fastq1="/path/to/input_read1.fastq"         # Path to the first paired-end FASTQ file
input_fastq2="/path/to/input_read2.fastq"         # Path to the second paired-end FASTQ file
sample_name="sample_identifier"                  # Sample identifier
known_sites_vcf="/path/to/known_sites.vcf"        # Known sites VCF file for BQSR (optional, if available)
output_directory="/path/to/output_directory"      # Directory for the output BAM file and reports
output_bam="${output_directory}/${sample_name}.bam"            # Output BAM file path
output_bqsr_report="${output_directory}/${sample_name}.BQSR-report.txt"  # Output BQSR report path

# Run the Parabricks fq2bam alignment pipeline
pbrun fq2bam \
    --ref ${genome_reference} \                # Reference genome (.fasta)
    --in-fq ${input_fastq1} ${input_fastq2} \  # Input paired-end FASTQ files
    --read-group-sm ${sample_name} \           # Sample name for read group
    --knownSites ${known_sites_vcf} \          # Known sites for BQSR
    --out-bam ${output_bam} \                  # Output BAM file
    --out-recal-file ${output_bqsr_report}     # Output Base Quality Score Recalibration (BQSR) report


```

Here is the command to submit job

```
chmod +x align.sh    # Makes the script executable
sbatch align.sh      # Submits the script to the SLURM queue
```

Use `squeue -u yourusername` to check job status.

## Additional Resources

### Datasets

In bioinformatics, it’s common to download databases or reference genomes from public websites. For example, performing sequence alignment requires downloading the appropriate reference genome. To simplify this process, we have pre-downloaded and managed several databases/datasets for users. These include:

- [Alphafold2](https://tuftsrt.github.io/guides/dev/bio/databases/doc/alphafold.html)

- [NCBI BLAST](https://tuftsrt.github.io/guides/dev/bio/databases/doc/ncbi.html)

- [Diamond](https://tuftsrt.github.io/guides/dev/bio/databases/doc/diamond.html)

- [iGenomes](https://tuftsrt.github.io/guides/dev/bio/databases/doc/igenomes.html)

- [Kraken2](https://tuftsrt.github.io/guides/dev/bio/databases/doc/kraken2.html)

- [Biobakery](https://tuftsrt.github.io/guides/dev/bio/databases/doc/biobakery.html)

- [geNomand](https://tuftsrt.github.io/guides/dev/bio/databases/doc/genomad_db.html)

- [Metaphlan](https://tuftsrt.github.io/guides/dev/bio/databases/doc/metaphlan.html)

### [New user guide](https://tuftsrt.github.io/guides/dev/bio/index.html)

In early 2025, we launched a new [RT Guides website](https://tuftsrt.github.io/guides/dev/index.html), offering comprehensive resources on a wide range of topics, including but not limited to HPC, data science, and and, **most importantly, [bioinformatics](https://tuftsrt.github.io/guides/dev/bio/index.html)**. We keep up with the latest trends and regularly update our materials to reflect new developments. We highly recommend bookmarking the website and referring to it whenever you encounter challenges. Your feedback is invaluable—let us know if you spot any errors or have suggestions.

For updates on bioinformatics education, software, and tools, consider [subscribing](https://elist.tufts.edu/sympa/info/best) to our e-list: [best@elist.tufts.edu](mailto:best@elist.tufts.edu).
