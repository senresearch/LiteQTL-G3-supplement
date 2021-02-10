This repository contains code that reproduces results from paper: [Speeding up eQTL scans in the BXD population using GPUs](https://www.biorxiv.org/content/10.1101/2020.06.22.153742v1.full.pdf)

- memory requirement: running the hippocampus data requires 81 GB of RAM. 
- Dependencies:
  - CUDA v10 
    - It is not necessary to have CUDA installed. There is a version that will run without access to CUDA. However, if you would like to see the timing comparison between CPU and GPU runs, Nvidia GPU is required, as well as CUDA library. Installation of CUDA library is not included in this README. Please refer to [Nvidia website](https://docs.nvidia.com/cuda/index.html) for guide to install CUDA library based on your operating system. 
  - Julia v1.4.0 
    - Julia packages: DelimitedFiles, Statistics, DataFrames, CSV, PyCall, PyPlot, LinearAlgebra, Random, Dates, CUDA, LiteQTL
  - R v4.0.2
    - R packages: tidyverse,qtl,mice,tictoc,qtl2, parallel  


How to run LiteQTL:  
- In your terminal, if you are at the the current repository of `LiteQTL-G3-supplement`, type `cd code`  
- If you have CUDA library installed on your machine, and would like to see the comparison between cpu runs and gpu runs, type `make all`. This command will 
  - install R and Julia dependency packages
  - download BXD dataset (including spleen and hippocampus phenotype data, and genotype data from genenetwork.org
  - clean phenotype data with R, tidyverse package, combine genotype and phenotype data to form a cross object, which r/qtl take as an input. 
  - run r/qtl to generate genotype probability and run genome scan
  - run julia genome scan with CPU and GPU, with Float32 and Float64 precision
  - compare the result agains r/qtl2 genome scan result. 
  - run a test of gemm, compare the timing of CPU and GPU, generate a figure to visualize the speed up gained from using GPU. 
- if you don't want to run gpu part of the code, type `make nogpu`, this command will include all steps in `make all` except installing julia's CUDA dependency package, and genome scan using Julia GPU, and gemm test. 


This repository is organized as follow: 
- code: This folder contains all necessary scripts to regenerate results from the paper, which includes a) downloading dataset used in the paper and put it in `data` folder, b) get dependency R and Julia packages, c) clean and prepare data for genome scan, d) run genome scan. A makefile is used to automate the process. 

- data: This folder contains all data used in the paper. Sub folder includes `raw` and `processed`. downloaded data will be placed in `raw`, and cleaned data will be placed in `processed`. Generated genome scan results will be placed in `results`.

- figures: generated figures will be placed here. 


[LiteQTL](https://github.com/ChelseaTrotter/LiteQTL.jl) is a julia package that uses GPU to accelerate eQTL genome scan using linear model.


Dataset:
- Genotype: BXD Genotype Database (GN Accession: GN600)
- Spleen: UTHSC Affy MoGene 1.0 ST Spleen (GN Accession : GN283)
- Hippocampus: UMUT Affy Hippocampus Exon (GN Accession: GN206)

