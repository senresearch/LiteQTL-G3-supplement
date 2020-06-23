This repository contains code that reproduces results from [LMGPU](https://www.biorxiv.org/content/10.1101/2020.06.22.153742v1.full.pdf)

[LMGPU](https://github.com/ChelseaTrotter/LMGPU.jl) is a julia package that uses GPU to accelerate eQTL genome scan using linear model. 
- It uses Julia 1.0.5
- Dependencies: 
  - [CUDA](https://developer.nvidia.com/cuda-92-download-archive), [CuArrays](https://github.com/JuliaGPU/CuArrays.jl), [CUDAnative](https://github.com/JuliaGPU/CUDAnative.jl), [CUDAdrv](https://github.com/JuliaGPU/CUDAdrv.jl)
  - Julia libraries: DelimitedFiles, LinearAlgebra
  
