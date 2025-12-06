# Organizing a Bioinformatics Project on the Tufts HPC

*For PhD students, rotation students, and beginners starting their first HPC project*             
Author: Shirley Li, xue.li37@tufts.edu             
Date: 2025-12-06         


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



## 6. Essential Best Practices 
- **Keep raw and processed data separate.** Treat raw data as read-only; archive it and move it to cold storage after the paper is published or the project is complete. For public datasets, document how you downloaded and preprocessed them, and delete raw public files once processing is finished—they can always be re-downloaded and should not occupy long-term storage.
- **Document software, tools, and package versions.** Record versions and parameters as you work, and write analysis notes in the style of a methods section rather than waiting until manuscript writing. Use the secure LLM tools available on the Tufts HPC to help organize and polish the wording.
- **Maintain reasonable folder sizes.** A typical manuscript folder (excluding raw and intermediate files) is small—often only a few GB. A full PhD project folder (excluding raw/intermediate data) should also remain within a manageable size. Check regularly; decluttering continuously is much easier than cleaning up everything at the end or when lab storage becomes tight.
- **Use consistent naming conventions.** Match analysis scripts with their SLURM submission files (e.g., `01_qc.R` ↔ `01_qc.sbatch`), and use logical, descriptive naming to keep workflows clear and easy to navigate.
- **Save intermediate processed objects.** Store Seurat `.rds` files or other processed objects to avoid re-running costly QC, loading, and preprocessing steps, especially for large-scale omics analyses where preprocessing can be time-consuming.
- **Track your code and documentation with Git/GitHub.** Version-control only code, notebooks, SLURM files, and documentation. Exclude all data via `.gitignore`, and sync regularly to maintain a clean, reproducible workflow history.
