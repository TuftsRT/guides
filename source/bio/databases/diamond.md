# Diamond: accelerated BLAST compatible local sequence aligner

DIAMOND is a sequence aligner for protein and translated DNA searches, designed for high performance analysis of big sequence data. The key features are:

- Pairwise alignment of proteins and translated DNA at 100x-10,000x speed of BLAST.
- Protein clustering of up to tens of billions of proteins
- Frameshift alignments for long read analysis.
- Low resource requirements and suitable for running on standard desktops or laptops.
- Various output formats, including BLAST pairwise, tabular and XML, as well as taxonomic classification.

Home page: https://github.com/bbuchfink/diamond

Tufts TTS Research Technology has built diamond database for NCBI non-redundant protein sequence (nr) database. You can use it to search your query sequences against nr using diamond.

- Update date: 04/05/2024
- Directory at HPC: `/cluster/tufts/biocontainers/datasets/diamond`
