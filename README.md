This repository contains code that reproduces results from paper: [Speeding up eQTL scans in the BXD population using GPUs](https://www.biorxiv.org/content/10.1101/2020.06.22.153742v1.full.pdf)

This repository is organized as follow:
- code: This folder contains all necessary scripts to regenerate results from the paper, which includes a) downloading dataset used in the paper and put it in `data` folder, b) get dependency R and Julia packages, c) clean and prepare data for genome scan, d) run genome scan. 

- data: This folder contains all data used in the paper. Sub folder includes `raw` and `processed`. downloaded data will be placed in `raw`, and cleaned data will be placed in `processed`. 

- figures: generated figures will be placed here. 


[LMGPU](https://github.com/ChelseaTrotter/LMGPU.jl) is a julia package that uses GPU to accelerate eQTL genome scan using linear model.
- It uses Julia 1.0.5
- Dependencies:
  - [CUDA](https://developer.nvidia.com/cuda-92-download-archive), [CuArrays](https://github.com/JuliaGPU/CuArrays.jl), [CUDAnative](https://github.com/JuliaGPU/CUDAnative.jl), [CUDAdrv](https://github.com/JuliaGPU/CUDAdrv.jl)
  - Julia libraries: DelimitedFiles, LinearAlgebra, BenchmarkTools,


Dataset:
you can run `make getdata` to download data. 
- Spleen (citation)
- Hippocampus (citation)

Data cleaning scripts:
- Dependencies: mice, r/qtl (citation), r/qtl2 (citation), tidyverse

Comparison:
- gemma (citation)
  - run script
- r/qtl2 (citation)
  - run script

Run gemm:
- run script
- Dependencies: LinearAlgebra, CuArrays, BenchmarkTools
