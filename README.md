# Nextflow Variant Calling Pipeline

## Overview
This repository contains a beginner-friendly and modular Nextflow DSL2 pipeline for variant calling from paired-end sequencing data.

The pipeline automatically performs reference genome indexing, quality control of raw reads, adapter and quality trimming, read alignment to a reference genome, and variant calling. It is designed for learning purposes, modular workflow understanding, and academic submissions.

---

## Pipeline Steps
1. **Reference Indexing (BWA Index)**  
   The reference genome FASTA file is indexed using BWA.  
   This step prepares the reference for efficient read alignment and is executed once per reference genome.

2. **FastQC – Quality Control**  
   Performs quality control analysis on raw paired-end FASTQ files to assess read quality.

3. **Read Trimming (fastp)**  
   Adapter removal and quality trimming of paired-end reads to improve alignment accuracy.

4. **Alignment (BWA-MEM + Samtools)**  
   Trimmed reads are aligned to the indexed reference genome, and sorted BAM files are generated.

5. **Variant Calling (VarScan)**  
   Variants are identified from aligned reads using mpileup-based variant calling, producing a VCF file.

---

## Directory Structure
nextflow_variant_pipeline/
├── data/
│ ├── samples/
│ │ ├── sample1_1.fastq
│ │ └── sample1_2.fastq
│ └── reference/
│ └── ref.fasta
├── modules/
│ ├── fastqc.nf
│ ├── trimming.nf
│ ├── alignment.nf
│ └── variant_calling.nf
├── main.nf
├── nextflow.config
└── README.md


---

## Requirements

The following software and system requirements are needed to run the pipeline:

- Linux or macOS (tested on WSL Ubuntu)
- Java (version 8 or higher)
- Nextflow (version 23.x or later)
- Basic command-line knowledge

## Tools Used

The pipeline integrates the following bioinformatics tools:

- **FastQC** – Quality control of raw sequencing reads
- **fastp** – Adapter removal and quality trimming
- **BWA** – Reference genome indexing and read alignment
- **Samtools** – BAM file sorting, indexing, and mpileup generation
- **VarScan** – Variant calling from mpileup files



## How to Run
From the project root directory:

```bash
nextflow run main.nf
Container and Resources
Executor: Local

CPUs: 2

Memory: 2 GB

Docker image: biocontainers/bcftools:v1.17-1-deb_cv1

Purpose of This Project
Learn Nextflow pipeline development

Understand modular workflow design

Practice basic variant calling concepts

Create a GitHub-ready academic project

Author
Harshvardhan Shetgar
MSc Bioinformatics Aspirant


---

If you want, next we can:
- Fix and clean `main.nf`
- Replace dummy variant calling with real `bcftools`
- Make the pipeline fully production-style (channels, params, logging)


Just tell me the next step and we will finish this properly.
