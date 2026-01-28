# Nextflow Variant Calling Pipeline

## Overview
This repository contains a simple and beginner-friendly variant calling pipeline built using Nextflow. The pipeline processes paired-end FASTQ files, performs basic quality control and preprocessing, aligns reads to a reference genome, and produces a VCF file.

This project is intended for learning purposes, modular workflow understanding, and academic submission.

---

## Pipeline Steps
1. FastQC – Quality control of raw sequencing reads  
2. Read Trimming – Adapter and quality trimming using Cutadapt  
3. Alignment – Read alignment to a reference genome using Minimap2  
4. Variant Calling – Generation of VCF output (simplified for learning)

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
- Nextflow
- Docker
- Linux or WSL environment

All bioinformatics tools are executed inside Docker containers to ensure reproducibility.

---

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