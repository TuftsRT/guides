---
tags: bioinformatics
---

# Parallel read alignment

Shirley Li, Bioinformatician, TTS Research Technology
xue.li37@tufts.edu

Date: 2024-11-01

This documentation outlines the process for running read alignment in parallel using the SLURM job scheduler. Specifically, it describes how to map and align Illumina reads to a reference genome by leveraging SLURM's job arrays and job submission capabilities. The approach involves using a wrapper script to read sample information from a text file and submit SLURM jobs for each sample.

## Prerequisites

Before proceeding, ensure you have the following:

- Access to a SLURM-managed high-performance computing (HPC) cluster.
- Installed modules for BWA and SAMtools (these tools have been installed in Tufts HPC as modules).
- A text file (samples.txt) containing the list of sample pairs.
- Read and write permissions for the data and output directories.

## Step-by-Step Guide

1. Prepare the Sample List
   Create a text file named `samples.txt`(or any other name) with each line containing the paired-end read files for a sample. The format should be:

```
sample1_R1.fastq sample1_R2.fastq
sample2_R1.fastq sample2_R2.fastq
...
sample20_R1.fastq sample20_R2.fastq
```

2. SLURM Script: `alignment_script.sh`
   Create a SLURM script that accepts read1 and read2 as arguments and performs the read alignment.

```bash
#!/bin/bash -l
#SBATCH -J parallel_alignment
#SBATCH -p batch
#SBATCH --n 4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=MyJob.%j.out
#SBATCH --error=MyJob.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=name@tufts.edu

# Load necessary modules
module load bwa/0.7.17
module load samtools/1.9

# Read the sample information from the arguments
READ1=$1
READ2=$2

# Define directories
DATA_DIR="/path/to/data"
REF_GENOME="/path/to/reference/genome.fa"
OUTPUT_DIR="/path/to/output"

# Define input files
READ1=${DATA_DIR}/${READ1}
READ2=${DATA_DIR}/${READ2}

# Define output files
OUTPUT_BAM=${OUTPUT_DIR}/sample.${READ1}.${READ2}.bam

# Align reads to the reference genome
# Adjust the parameters according to the best practices in your field
bwa mem -t 4 ${REF_GENOME} ${READ1} ${READ2} | samtools view -bS - > ${OUTPUT_BAM}

# Optionally, you can sort and index the BAM file
samtools sort -o ${OUTPUT_DIR}/sample.${READ1}.${READ2}_sorted.bam ${OUTPUT_BAM}
samtools index ${OUTPUT_DIR}/sample.${READ1}.${READ2}_sorted.bam

echo "Alignment for sample ${SLURM_ARRAY_TASK_ID} completed"

```

3. Wrapper Script: submit_jobs.sh
   Create a wrapper script to read the sample list and submit jobs to SLURM.

```bash
#!/bin/bash -l

# Define variables
SAMPLES_FILE="samples.txt"
SCRIPT="alignment_script.sh"

# Submit the job for each sample
# The while loop will read each line to find READ1 and READ2
while read READ1 READ2
do
  sbatch ${SCRIPT} ${READ1} ${READ2}
done < $SAMPLES_FILE

```

4. Execute the Scripts
   4.1 Ensure both scripts are executable:

```
chmod +x alignment_script.sh
chmod +x submit_jobs.sh
```

4.2 Run the wrapper script to submit the jobs:

```
./submit_jobs.sh
```

## Conclusion

By following this documentation, you can efficiently run read alignment for multiple samples in parallel on an HPC cluster using SLURM. This approach optimizes resource usage and reduces the overall processing time for large datasets.
