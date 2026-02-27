# Running Bioinformatics Jobs on Tufts HPC

2026-02-24           

Shirley Li: [xue.li37@tufts.edu ](mailto:xue.li37@tufts.edu)         



## Workshop Overview

This hands-on workshop introduces how to run bioinformatics analyses on the Tufts HPC using SLURM. Using an RNA-seq example, participants will learn how to structure projects, submit jobs, and scale workflows efficiently.

**You will learn how to:**

- Organize a reproducible HPC project directory

- Write and submit SLURM batch scripts

- Allocate CPU, memory, and runtime appropriately

- Monitor logs and troubleshoot jobs

- Use job dependencies and SLURM job arrays for scalable analysis

  

## 1. Set Up Your Working Directory

All work for this project will be done under:

```
/cluster/tufts/workshop/utln/
```

If you have a dedicated lab storage space, you may use that instead.



Create a directory for your project

```
mkdir /cluster/tufts/workshop/utln/myproject/    
# creating a folder to hold all files related to this analysis
```

**Always create a dedicated project directory. Never mix analyses together.**



## 2. Create a Basic Project Structure

Navigate into your project directory and create the following folders:

```
cd /cluster/tufts/workshop/utln/myproject/
# You are now inside your project folder

mkdir raw_data results scripts
touch README.md 
# Document what this project does

mkdir scripts/logs
# Create logs folder under scripts
# SLURM will generate: .out files (standard output) and .err files (error messages)
# Keeping logs in one place makes debugging much easier.
```



Your project structure should look like this:

```
myproject/
├── raw_data/        # input files (never modify these)
├── results/         # output files from analysis
├── scripts/         # your SLURM scripts
│   └── logs/
└── README.md
```

This structure keeps:

- Raw input data separate
- Results organized
- Scripts and logs clearly managed



## 3. Prepare Input Fastq Files

The example FASTQ files are located at:

```
/cluster/tufts/workshop/public/2026spring/nfcore/fastq/*
```

Instead of copying large files, create symbolic links in your `raw_data` directory:

```
cd /cluster/tufts/workshop/utln/myproject/raw_data/
# move into raw_data folder

ln -s /cluster/tufts/workshop/public/2026spring/nfcore/fastq/* ./
# Create symbolic links
# This creates shortcuts in your raw_data/ folder that point to the original files.
```



**Why We Don’t Copy the Files**

FASTQ files are large. Copying them would:

- Waste storage
- Slow down the system
- Create unnecessary duplicates

Instead, we use symbolic links.



This allows you to:

- Access the files locally in your project

- Avoid duplicating large datasets

- Keep the original data unchanged

  

## 4. Running FastQC

We will now write our first batch job to run FastQC.

`fastqc.sh`

```
#!/bin/bash
#SBATCH --job-name=fastqc_rnaseq
#SBATCH --partition=preempt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/fastqc_%j.out
#SBATCH --error=logs/fastqc_%j.err

echo "Job started on $(date)"
echo "Running on node: $(hostname)"

# Load FastQC module (adjust if needed)
module load fastqc/0.11.9  


DIR=/cluster/tufts/workshop/utln/myproject/
# change this line

mkdir -p $DIR/results/fastqc

# Run FastQC on all FASTQ files
fastqc $DIR/raw_data/*.fastq.gz \
       --outdir $DIR/results/fastqc \
       --threads 8

echo "Job finished on $(date)"
```

**What This Job Does**

- Reads FASTQ files from `raw_data/`
- Writes reports to `results/fastqc/`
- Saves logs in `scripts/logs/`
- Uses 8 CPU threads
- Runs all samples in one submission

Submit:

```
sbatch fastqc.sh
```



## 5. Alignment with STAR

`star.sh`

```
#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --cpus-per-task=8
#SBATCH --partition=preempt
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/star_%j.out
#SBATCH --error=logs/star_%j.err

module load star/2.7.11b     

STARINDEX=/cluster/tufts/workshop/public/2026spring/star_index/
DIR=/cluster/tufts/workshop/utln/myproject/

mkdir -p $DIR/results/star/

for fq in $DIR/raw_data/*_1_*.fastq.gz
do
    sample=$(basename $fq _1_sub.fastq.gz)
    
    STAR \
      --runThreadN 8 \
      --genomeDir $STARINDEX \
      --readFilesIn $fq \
      --readFilesCommand zcat \
      --outFileNamePrefix $DIR/results/star/${sample}_
      
    echo $fq DONE
done
```

**What This Job Does**

- Aligns reads to the reference genome
- Processes all samples in a loop
- Writes BAM output to `results/star/`
- Uses 8 CPU threads

Submit:

```
sbatch star.sh
```



## 6. Post-Processing the BAM File (Sorting)

Many downstream tools require a **sorted BAM file**. We will sort the BAM file using `samtools`.

`sort.sh`

```
#!/bin/bash
#SBATCH --job-name=sort_bam
#SBATCH --cpus-per-task=4
#SBATCH --partition=preempt
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/sort_%j.out
#SBATCH --error=logs/sort_%j.err

module load samtools/1.21 

mkdir -p $DIR/results/sorted_bam/

for bam in $DIR/results/star/*_Aligned.out.bam
do
    sample=$(basename $bam _Aligned.out.bam)

    samtools sort \
        -@ 4 \
        -o $DIR/results/sorted_bam/${sample}_sorted.bam \
        $bam
done

```

**What This Job Does**

- Takes STAR output BAM files
- Sorts each BAM file
- Writes sorted files to `results/sorted_bam/`
- Uses 4 CPU threads

Submit:

```
sbatch sort.sh
```



## 7. Wrapper Script to Chain Jobs

This script chains STAR and sorting using SLURM dependency. FastQC runs independently.

`run_pipeline.sh`

```
#!/bin/bash

# Submit FASTQC independently
sbatch fastqc.sh

# Submit STAR alignment
jid1=$(sbatch star.sh | awk '{print $4}')
echo "STAR job submitted with Job ID: $jid1"

# Submit sorting job after STAR completes successfully
jid2=$(sbatch --dependency=afterok:$jid1 sort.sh | awk '{print $4}')
echo "Sorting job submitted with Job ID: $jid2"

```



Make Script Executable and run:

```
chmod +x run_pipeline.sh
./run_pipeline.sh
```



**What This Does**

- FASTQC runs independently
- STAR runs
- Sorting runs only if STAR finishes successfully
- If STAR fails, sorting will not start



------



## Advanced usage: SLURM Job Array

In the previous `star.sh` script, we used a `for` loop to process all samples inside a single SLURM job.

That approach works, but it runs samples **sequentially** — one after another — within the same job allocation.

If you have multiple independent samples, a more scalable approach is to use a **SLURM job array**.



### Why Use a Job Array?

A job array allows SLURM to:

- Run multiple jobs of the same script
- Process different input files independently
- Execute jobs in parallel across different compute nodes

Instead of one job looping through 6 samples, SLURM launches 6 separate jobs automatically.

Each job:

- Processes one sample
- Has its own log file
- Uses its own allocated resources

This improves efficiency and scalability.

`star_array.sh`

```
#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --partition=preempt
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-6
#SBATCH --output=logs/star_%A_%a.out
#SBATCH --error=logs/star_%A_%a.err

module load star/2.7.11b 

STARINDEX=/cluster/tufts/workshop/public/2026spring/star_index/
DIR=/cluster/tufts/workshop/utln/myproject/

mkdir -p $DIR/results/star_arrayjob/

#============================
# --array=1-6: where 6 = number of samples.
# If you don’t know how many samples there are, you can count them first:
# ls $DIR/raw_data/*_1_*.fastq.gz | wc -l
#============================

files=($DIR/raw_data/*_1_*.fastq.gz)
# This collects all matching FASTQ files into a bash array.

fq=${files[$SLURM_ARRAY_TASK_ID-1]}
# SLURM_ARRAY_TASK_ID starts at 1.
# Bash arrays start at 0.

sample=$(basename $fq _1_sub.fastq.gz)

STAR \
  --runThreadN 8 \
  --genomeDir $STARINDEX \
  --readFilesIn $fq \
  --readFilesCommand zcat \
  --outFileNamePrefix $DIR/results/star/${sample}_
```



