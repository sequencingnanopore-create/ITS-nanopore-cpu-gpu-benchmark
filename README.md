
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.19325905)

## Overview
This repository contains all bioinformatics scripts used in the study:
> 'Machine Learning–Enhanced Nanopore ITS Analysis: Evaluating CPU–GPU Pipelines for High-Accuracy Fungal Taxonomic Resolution'
> Albuja et al. (2025), Frontiers in Bioinformatics.

## System Requirements
- CPU workflow: Google Colaboratory (free tier) or Linux with Python 3.8+
- GPU workflow: HPC with NVIDIA GPU, CUDA ≥ 11.0, Dorado v0.9.6

## Installation
```bash
git clone https://github.com/[user]/ITS-nanopore-cpu-gpu-benchmark.git
cd ITS-nanopore-cpu-gpu-benchmark
pip install -r requirements_cpu.txt
```

## Data Availability
Raw sequencing data are deposited in NCBI SRA under BioProject PRJNA[XXXXXX].

## HPC Basecalling Commands (Dorado SUP)
The GPU basecalling was performed on CEDIA HPC (NVIDIA A100-SXM4-40GB).
Commands are documented in docs/HPC_basecalling_commands.md

## Citation
If you use this code, please cite: [reference]
