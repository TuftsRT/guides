Date: Feb 6, 2024\
Author: Shirley Li, xue.li37@tufts.edu\
Class ID: Sp24-IDGH-1001-1-Bioinformatics\
Canvas link: https://canvas.tufts.edu/courses/55751

# Introduction to Metagenomics - Session 2

## Learning Objective

1. NCBI Database Proficiency: Develop skills to efficiently locate and interpret data on the [NCBI database](https://www.ncbi.nlm.nih.gov/), including navigating to specific BioProject and SRA experiment pages.
1. Data Retrieval from Published Papers: Gain the ability to identify and extract relevant raw data and metadata from published scientific papers.
1. Metagenomic Sequencing Platforms: Learn about different sequencing platforms by analyzing their data characteristics, specifically focusing on [Illumina](https://www.sciencedirect.com/science/article/pii/S0198885921000628) and [Nanopore technologies](https://www.nature.com/articles/nbt.3423).
1. Taxonomic Analysis: Acquire practical experience in assigning taxonomic labels to sequencing reads using Kraken2 on [Tufts Galaxy](https://galaxy.cluster.tufts.edu/), and in converting and visualizing these labels with Krona.
1. Data Comparison and Interpretation: Enhance skills in comparing visualized data with NCBI SRA information and drawing conclusions about sample composition.

## Exercise 1: NCBI Database Navigation

### Objective:

In this section, you will learn to navigate the [NCBI database](https://www.ncbi.nlm.nih.gov/bioproject). Known for its comprehensive nature, the database hosts a variety of sections, each tailored to specific types of datasets. Your focus will be on mastering the ability to identify and navigate through BioProject pages, SRA experiment pages, and SRA Runs linked to published research papers. This proficiency is crucial in the field of bioinformatics, enabling us to effectively leverage previous research as a foundation for new discoveries and advancements in the study of biological data.

<img src="https://github.com/shirleyxueli41/Tufts_workshops/assets/88347911/32383594-b895-4f0b-9107-882313c69304" width="900" height="300">

<em>A screenshot of NCBI website</em>

### Instructions:

Your task involves exploring the NCBI database to gather specific information:

1. Review the research paper ["Evaluation of full-length nanopore 16S sequencing for detection of pathogens in microbial keratitis"](https://peerj.com/articles/10778/) . Your goal is to identify the associated BioProject ID within the paper.
   - Hint: BioProject ID is [PRJEB37709](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB37709/) in the "Data Availability" section of the paper.
1. Answer the following queries regarding the BioProject, and keep a record of your findings:\
   a. What is the URL for this specific BioProject?\
   b. Total number of biosamples included in this BioProject.\
   c. Total number of SRA experiments associated with this BioProject.\
   d. Determine the sequencing platform used for SRA run [ERR4836970](https://www.ncbi.nlm.nih.gov/sra/?term=ERR4836970).
1. Apply the same procedure as step 2 for a second research paper: ["Benchmarking second and third-generation sequencing platforms for microbial metagenomics"](https://www.nature.com/articles/s41597-022-01762-z).
   - Hint: BioProject ID is in the "Data Records" section of the paper.

## Exercise 2: NCBI Database Exploration

### Objective:

Engage in a hands-on exercise to explore the NCBI database using specific SRA run IDs. Your task will involve navigating various sections of the database and applying your understanding of sequencing technologies to hypothetical research scenarios.

### Instructions:

Utilize the SRA run ID to search the NCBI website. Explore the corresponding SRA, BioSample, and BioProject sections related to this SRA run.

1. Assignment Completion: Choose one SRA run and document your findings in the provided Google Spreadsheet: [Exercise Spreadsheet](https://docs.google.com/spreadsheets/d/1s_NVSPLABQtTmB-EXwY4gmNOzGJxr2RE7WVnbFz6djw/edit#gid=927665114).

<img src="https://github.com/shirleyxueli41/Tufts_workshops/assets/88347911/de1556b2-f038-451a-b4dc-a4f3d324aeae" width="900" height="300">

<em>A screenshot of the spreadsheet</em>

2. Questions for Analysis\
   a. **Mars Soil Sample Analysis**: If you obtained a soil sample from Mars for identifying microorganisms and assembling their genomes, which sequencing technology would be optimal? Consider factors like the detection of novel organisms and the precision required for genome assembly. Discuss your choice, focusing on read length, accuracy, and cost implications.\
   b. **Gut Microbiome Study**: In researching the impact of dietary changes on the gut microbiome, what type of sample would you collect, and which sequencing technology would be most suitable? Provide your rationale for this choice.

Additional Resources: Hints for these questions can be found [here](https://github.com/shirleyxueli41/Tufts_workshops/blob/main/IGDH-1001_2024Feb/Exercise%202_hints.pdf).

## Exercise 3 Taxonomy assignment and interpretation.

### Objective:

Use Kraken2 for taxonomy assignment and visualize the results with a Krona plot. Interpret and present the result.

### Instructions:

> [!NOTE]
> The tools we will use for this analysis are:
>
> - Download and Extract Reads in FASTA/Q
> - Kraken2
> - filter
> - sort

1. Log in to your [Galaxy account](https://galaxy.cluster.tufts.edu/).

1. Name the history as "Session 2 Metagenomics-ERR12302112" by double clicking the "Unnamed history".

1. Now let's start the analysis:

   1. Under tools on the far left of the page, search for **Download and Extract Reads in FASTA/Q** format from NCBI SRA, run the tool with the following parameters:\
      \- **Accession:** ERR12302112\
      \- Click **Execute**\
      ![Screenshot 2024-01-25 at 16 09 58](https://github.com/shirleyxueli41/Tufts_workshops/assets/88347911/9cafa3e2-78ad-4040-9e7f-a3baae378512)

   1. **Kraken2** assign taxonomic labels to sequencing reads with the following parameters:\
      \- **Single or paired end**: Single\
      \- **Input Sequences**: the output from last step. Ex: 1.ERR12302112 (fastq-dump)\
      \- Click Create Report, then set **Print a report with aggregate counts/clade to file** to Yes\
      \- **Select a Kraken2 database**: Minikraken2 v2\
      ![Screenshot 2024-01-25 at 16 17 44](https://github.com/shirleyxueli41/Tufts_workshops/assets/88347911/26d7e312-a8a7-49c5-b3d5-39e9ce5893d9)\
      _Note this step will create two output files_\
      <img width="296" alt="Screenshot 2024-01-25 at 16 21 42" src="https://github.com/shirleyxueli41/Tufts_workshops/assets/88347911/492691e2-efe6-49da-9368-ab33c1312d21">

   1. **Filter** data on any column using simple expressions with the following parameters:\
      \- **Filter**: the _report_ output from last step. Ex: Report: Kraken2 on data 1\
      \- **With following condition**: c4=="S"\
      _This will keep the rows whose fourth column has a character S, S stands for species_

   1. **Sort** data in ascending or descending order with the following parameters:\
      \- **Sort Dataset**: the output file from filter. Ex: Filter on data 2\
      \- **with flavor**: Numerical sort\
      \- **everything in**: Descending order\
      _Take a look at the output file, the first few lines should be like this:_

      <img width="743" alt="Screenshot 2024-01-25 at 16 30 26" src="https://github.com/shirleyxueli41/Tufts_workshops/assets/88347911/d284b4de-d6db-43fb-b14b-5c2452823902">

### In-class assignment:

Divide into teams (either two or three teams). Each team should select one SRA run from the provided [google spreadsheet](https://docs.google.com/spreadsheets/d/1s_NVSPLABQtTmB-EXwY4gmNOzGJxr2RE7WVnbFz6djw/edit#gid=941310028). Then, replicate the previously outlined steps to identify the top three prevalent species. Research one or two of these species using Google, and compare your findings with the samples to check for coherence. Each team will be given five minutes to showcase their findings. An example report can be found [here](https://github.com/shirleyxueli41/Tufts_workshops/blob/main/IGDH-1001_2024Feb/Exercise%203_Example%20Report.pdf).

> [!WARNING]\
> _Warning: Ensure you generate a fresh history and assign a distinct name for the analysis._\
> Click the "+" button on the top right to create new history session.

## Exercise 4 Taxonomy visualization.

### Objectives:

The exercise aims to utilize Krona for creating interactive visualizations of taxonomic data, highlighting the tool's effectiveness in representing complex hierarchical structures. It also involves a comparison with NCBI Kroa, assessing differences in visualization techniques and data representation.

### Instructions:

1. Switch back to the session "Session 2 Metagenomics-ERR12302112".
1. **Krakentools: Convert kraken report file** to krona text file with the following parameters:
   - **Kraken report file**: The _report_ output from _Kraken2_. Ex: Report: Kraken2 on data1\
     _This will generate an output called "Krakentools: Convert kraken report file on data 2"_
1. **Visualize with Krona** Visualize any hierarchical data with the following parameters:
   - **Select input file:** Krakentools: Convert kraken report file on data 2.\
     _This will generate an output called "Krona on data 5: HTML"_\
     <img width="885" alt="Screenshot 2024-01-25 at 16 42 45" src="https://github.com/shirleyxueli41/Tufts_workshops/assets/88347911/8e1bd1c1-5708-47f7-9e0b-2ba4ac8f7559">
1. Compare the Krona plot with it on NCBI SRA. Link is [here](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=ERR12302112&display=analysis).
   - Click Show Krona View\
     _NCBI uses Sequence Taxonomic Analysis Tool ([STAT](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02490-0)), a scalable k-mer-based tool for fast assessment of taxonomic diversity intrinsic to submissions, independent of metadata._

# Reference

- https://bisonnet.bucknell.edu/files/2021/05/Kraken2-Help-Sheet.pdf
- https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00900-2
- https://jddtonline.info/index.php/jddt/article/view/5433
- https://www.sciencedirect.com/science/article/pii/S094450132200194X?via%3Dihub
- https://benlangmead.github.io/aws-indexes/k2
