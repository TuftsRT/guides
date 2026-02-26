# How to Organize a Bioinformatics Project on the Tufts HPC
*A Practical Guide for Reproducible and Scalable Research*
           
Author: Shirley Li, xue.li37@tufts.edu             
Date: 2026-02-26        

If you are running RNA-seq, single-cell, spatial transcriptomics, or other omics analyses on the Tufts HPC cluster, organizing your project correctly from the beginning will save time, reduce errors, and make your work reproducible.

This guide explains:
- Where to work on the Tufts HPC
- A recommended project folder structure
- What belongs in each directory
- How to use Git/GitHub correctly
- Best practices for long-term project management
- This structure works for SLURM-based workflows, R/Python pipelines, and multi-analysis projects.

## 1. Where to Work on Tufts HPC

Each user works under:

```
/cluster/tufts/<labname>/<utln>/
```

Your lab should maintain shared resources under:

```
/cluster/tufts/<labname>/shared/
```

Use the shared folder for reference genomes, indexes, annotation files, or datasets reused across the lab.



## 2. Recommended Project Structure 

Use this structure **inside your project**:

```
project_name/
    data_raw/
    data_processed/
    scripts_global/
    analysis1/
        scripts/
        results/
        logs/
        slurm/
        notebooks/
    analysis2/
        scripts/
        results/
        logs/
        slurm/
        notebooks/
    envs/
    docs/
    README.md
```

**Key idea:**
 Reusable code lives in `scripts_global/`; analysis-specific code stays inside each `analysis*/scripts/`.



## 3. What Goes in Each Folder 

### **data_raw/**

Untouched input files (FASTQ, BAM, MERFISH CSV, Visium outputs).
 Should never be modified.

### **data_processed/**

Filtered matrices, Seurat objects (`.rds`), processed tables.
 Time-consuming results that you want to reuse.

### **scripts_global/**

Shared utilities, helper functions, custom plotting, reusable QC code.

### **analysis\*/scripts/**

Analysis-specific R/Python scripts, e.g.:

```
01_qc.R
02_clustering.R
03_spatial.R
```

### **analysis\*/notebooks/**

Exploratory notebooks (Rmd/Jupyter).
 Not part of the formal pipeline.

### **analysis\*/slurm/**

Job submission files (`.sbatch`).
 Match names to scripts, e.g.:

```
scripts/01_qc.R  
slurm/01_qc.sbatch
```

### **analysis\*/results/**

Figures, tables, intermediate outputs.

### **analysis\*/logs/**

SLURM `.out`/`.err` logs and analysis logs.

### **envs/**

Conda environments, `renv` folder, Apptainer definitions.

### **docs/**

Notes, meeting summaries, methods drafts, descriptions of workflows.



## 4. Version Control: Use Git + GitHub

Track **only code and documents**, not large data.

Put these under Git:

```
scripts_global/
analysis*/scripts/
analysis*/slurm/
analysis*/notebooks/
docs/
README.md
```

Add to `.gitignore`:

```
data_raw/
data_processed/
results/
*.rds
*.h5
*.fastq*
*.bam
```

Create a GitHub repo for tracking scripts + documentation.


## 5. Essential Best Practices for HPC Projects
### 5.1. Separate Raw and Processed Data
Raw data:
- Never modified
- Archived after publication
- Moved to cold storage if appropriate
Public raw datasets:
- Document download commands
- Delete after processing (can be re-downloaded)

### 5.2. Document Software Versions
Record:
- Tool versions
- Parameters
- Package versions
Write analysis notes in a methods-style format while working.

### 5.3. Keep Folder Sizes Reasonable
A typical manuscript project (excluding raw/intermediate data) should only be a few GB.
Regular cleanup is much easier than emergency cleanup.

### 5.4. Use Consistent Naming Conventions
Match:
`01_qc.R  ↔  01_qc.sbatch`

Use logical numbering.
Avoid ambiguous names like:
`final_new_v2_fixed.R`

### 5.5. Save Intermediate Objects
For large omics analyses:
- Save Seurat .rds
- Save processed matrices
Preprocessing can take hours. Avoid repeating it.

### 5.6. Track Everything with Git

Track:
- Code
- SLURM scripts
- Notebooks
- Documentation
- Exclude all data.
- Commit frequently.

## Frequently Asked Questions
### Q: Should each paper have its own project folder?
Yes.
Each manuscript should have its own top-level project directory. This prevents cross-contamination of scripts, results, and documentation between studies.

### Q: Should I organize by data type or biological question?
Organize by biological question or manuscript goal.

For example:
```
analysis1/  → differential expression
analysis2/  → spatial analysis
analysis3/  → integration
```

Avoid mixing unrelated analyses in the same folder.

### Q: Where should I run interactive work on Tufts HPC?

Use Open OnDemand for:
- RStudio
- Jupyter
- Small-scale exploratory work

Use SLURM batch jobs for:
- Large RNA-seq runs
- STAR alignment
- nf-core workflows
- High-memory jobs

### Q: Should I run everything through SLURM?

If a job:
- Runs > 5–10 minutes
- Uses multiple cores
- Requires significant memory
It should be submitted via SLURM.

Interactive sessions are for exploration only.

### Q: Where should I store large reference genomes?

Store them in:
`/cluster/tufts/<labname>/shared/`

Never duplicate reference genomes across personal project folders.

### Q: How should I name my project folders?

Use descriptive and date-aware names:

```
2026_skin_merfish_project/
2026_rnaseq_inflammation/
```

Avoid vague names like:
```
new_project/
final_analysis/
test2/
```
### Q: Should notebooks replace scripts?
No.
Notebooks are for:
- Exploration
- Visualization
- Demonstration

Scripts are for:
- Reproducible pipelines
- SLURM jobs
- Formal analysis

A clean workflow uses both.

### Q: When should I clean up intermediate files?

After:
- You save key processed objects
- The pipeline has completed successfully
- The manuscript is accepted
- Large intermediate files (temporary BAMs, intermediate matrices) should not live forever.

### Q: How large should a project folder be?

As a rule of thumb:
Code + figures + docs → usually only a few GB

Large raw/intermediate data should not dominate your active project space

If a folder grows unexpectedly large, check:
`du -sh *`

### Q: Should I store environments inside the project?
Yes.
Keep environment definitions in:
`envs/`

This ensures:
- Reproducibility
- Easier collaboration
- Future reruns of analysis

### Q: How do I make my project reproducible for future lab members?

Keep clear folder structure
- Maintain README.md
- Track code with Git
- Record software versions
- Avoid undocumented manual steps

A future lab member should be able to understand your project in one hour.

### Q: What is the most common mistake in HPC project organization?

Mixing everything in one folder.
Typical bad pattern:
```
scripts/
data/
results/
more_scripts/
new_results/
final/
```
Without structure, scaling becomes impossible.

### Q: Should I compress FASTQ or BAM files?

FASTQ should remain .fastq.gz.
Do not unzip them.

BAM files can be archived if no longer needed, but do not delete primary alignment results before publication unless archived properly.
