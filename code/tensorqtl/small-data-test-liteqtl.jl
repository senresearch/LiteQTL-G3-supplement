using LiteQTL 
using CSV 
using DelimitedFiles
using DataFrames
using StatsBase
using CUDA
using Dates
using Base.Threads

using BenchmarkTools

# include("new-functions-liteqtl.jl")


pheno_file = joinpath(@__DIR__, "..", "..", "data", "tensorqtldata", "cleanpheno.csv")
geno_file = joinpath(@__DIR__, "..", "..", "data", "tensorqtldata", "chr9.csv")

pheno = LiteQTL.get_pheno_data(pheno_file, Float64, transposed=true)
# 10.865656 seconds (47.21 M allocations: 1.625 GiB, 12.64% gc time)
geno = CSV.read(geno_file, DataFrame)

datatype = Float32
# Only taking the first 5000 to get the results right. 
genomat = convert(Matrix{datatype}, geno[1:30000,3:end])|> transpose |> collect;
G = genomat


Y = convert(Matrix{datatype}, pheno)

LiteQTL.set_blas_threads(16);
n = size(Y,1)
m = size(Y,2)
p = size(G,2)

BenchmarkTools.DEFAULT_PARAMETERS.samples = 10 

println("******* Indivuduals n: $n, Traits m: $m, Markers p: $p ****************");
export_matrix = true
cputime = @benchmark lodc = scan(Y, G, export_matrix=export_matrix, usegpu=false, lod_or_pval="lod")
CUDA.@elapsed lodg = scan(Y, G, usegpu=true);
gputime = @benchmark CUDA.@sync lodg = scan(Y, G, usegpu=true, export_matrix=false);

open("time_compare_table.txt", "a") do io
    write(io, "$(now()) \t n:$n,m:$m,p:$p LiteQTL: \t nthreads: $(nthreads()) \t cpu(s):$(cputime.times[#=median=#2]) \t gpu(s)$(gputime.times[#=median=#2])\n")
end   

# writedlm("julia-scan-pval-result-cpu.csv", lodc, ',')


# covarfile = "../../data/tensorqtldata/covariates.csv"
# covar = CSV.read(covarfile, DataFrame)
# X = convert(Matrix{datatype}, covar[:, 2:end])|> collect


# julia> cputime
# BenchmarkTools.Trial:
#   memory estimate:  2.22 GiB
#   allocs estimate:  495
#   --------------
#   minimum time:     7.576 s (0.09% GC)
#   median time:      7.576 s (0.09% GC)
#   mean time:        7.576 s (0.09% GC)
#   maximum time:     7.576 s (0.09% GC)
#   --------------
#   samples:          1
#   evals/sample:     1

# julia> gputime
# BenchmarkTools.Trial:
#   memory estimate:  34.00 MiB
#   allocs estimate:  794
#   --------------
#   minimum time:     132.839 ms (0.00% GC)                                       
#   median time:      147.226 ms (2.21% GC)                                       
#   mean time:        148.532 ms (3.40% GC)                                       
#   maximum time:     168.271 ms (5.60% GC)                                       
#   --------------
#   samples:          10
#   evals/sample:     1