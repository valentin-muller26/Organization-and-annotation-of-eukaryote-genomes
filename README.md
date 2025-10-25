# Organization-and-annotation-of-eukaryote-genomes

## Description and goal of the project

This repositery contains scripts for the annotation of an assembly of *Arabidopsis thaliana* accession Lu-1. This assembly was done during the genome and transcriptome assembly course. 

## Datasets

**Assembly Data:**
The assembly for the genomic data come from the assesion Lu-1 and was sequence using Pacbio HiFi reads and assembled using Hifiasm version 0.25.0.

**Transcriptomic Data:**
The assembly for the transcriptomic data come from the assesion Sha and was sequence using Illumina paired-end short reads and assembled using Trinity version 2.15.1


The sequencing data can be found in the following link :  
- Genome data: https://www.nature.com/articles/s41588-024-01715-9
- Transcriptome data: https://www.nature.com/articles/s41467-020-14779-y

## TE Annotation Pipeline

### 1.a TE Annotation using EDTA `01a_run_EDTA.sh`

In this first step we used EDTA (Extensive de novo TE Annotator) version 2.2 to annotate the TEs in our assembly using the script `01a_run_EDTA.sh`. This annotation resulted in a folder containing TE library for the whole genome and classification of the elements into superfamilies.

More information about EDTA can be found [here](https://github.com/oushujun/EDTA)

### 1.b-d Analysing the Full-length LTR Retrotransposons 

Once EDTA has generated the annotation of the retrotransposons, we perform an analysis of the intact full-length LTR retrotransposons.

This analysis was performed in three steps:

**Step 1: Extract LTR percent identity `01b_extract_LTR_identity.sh`**

This script parses the main information of intact full-length LTR retrotransposons and provides some summary statistics.

**Step 2: Refine LTR-RT classification with TEsorter `01c_Run_TEsorter_FullLTR_classification.sh`**

This step uses TEsorter version 1.3.0 to refine the classification of the intact LTRs into clades.

**Step 3: Visualize LTR-RT dynamics `01d-full_length_LTRs_identity.R`**

This step uses R version 4.5.0 and R package `circlize` to generate circular plots showing the distribution and density of the Copia and Gypsy clades.

This analysis is done locally and the files `genomic.fna.mod.LTR.intact.raw.gff3` and `genomic.fna.mod.LTR.raw.fa.rexdb-plant.cls.tsv` need to be downloaded from the cluster and placed in the same directory as the `01d-full_length_LTRs_identity.R` script.


## List of the tools used

| Tool | Version | 
|------|---------|
| EDTA | 2.2 | 
| TEsorter | 1.3.0 |
| SeqKit | 2.0+ | 
| parseRM.pl | Latest | 
| circlize | R package |
| MAKER | 3.01.03 | 
| OpenMPI | 4.1.1 |
| BioPerl | 1.7.8 |
|R | 4.5.0|


## Author

Valentin Müller