using LiteQTL 
using CSV 
using DelimitedFiles
using DataFrames
using StatsBase
using CUDA
using Dates
using Base.Threads
using Statistics

using BenchmarkTools

# include("new-functions-liteqtl.jl")

timing_file = "/home/xiaoqihu/git/LiteQTL-G3-supplement/code/tensorqtl/liteqtl_timing_report.txt"
pheno_file = joinpath(@__DIR__, "..", "..", "data", "tensorqtldata", "cleanpheno.csv")
geno_file = joinpath(@__DIR__, "..", "..", "data", "tensorqtldata", "chr9.csv")

pheno = LiteQTL.get_pheno_data(pheno_file, Float64, transposed=true)
# 10.865656 seconds (47.21 M allocations: 1.625 GiB, 12.64% gc time)
geno = CSV.read(geno_file, DataFrame)

datatype = Float32
# Only taking the first 20000 to get the results right. 
genomat = convert(Matrix{datatype}, geno[1:20000,3:end])|> transpose |> collect;
G = genomat
Y = convert(Matrix{datatype}, pheno)

LiteQTL.set_blas_threads(Threads.nthreads());
n = size(Y,1)
m = size(Y,2)
p = size(G,2)

BenchmarkTools.DEFAULT_PARAMETERS.samples = 10 
println("******* Indivuduals n: $n, Traits m: $m, Markers p: $p ****************");

open(timing_file, "w") do io
    write(io, "time,device,data_transfer_time,compute_time,result_reorg_time,pval_calc_time,elapsed_total\n")
end 


############################# Full Matrix Case ###################################
cputime = 0 
gputime = 0
for i in 1:10
    export_matrix = true
    global cputime = @elapsed global lodcfull = scan(Y, G, maf_threshold=0.00, export_matrix=export_matrix, usegpu=false, lod_or_pval="lod", timing_file=timing_file)
    global gputime = CUDA.@elapsed global lodgfull = scan(Y, G,  maf_threshold=0.00, usegpu=true, export_matrix=export_matrix, timing_file=timing_file);
    # gputime = CUDA.@elapsed lodg = scan(Y, G, maf_threshold=0.0, usegpu=true, export_matrix=false);
end
@debug begin 
    "Results of compraring CPU full lod score to GPU full lod score matrix: $(sum(isapprox.(lodcfull, lodgfull, atol=0.001)) == prod(size(lodcfull)))"
end

open(timing_file, "a") do io
    write(io, "#$(now()),n:$n m:$m p:$p,LiteQTL Full Matrix,nthreads: $(nthreads()),cpu(s):$(cputime),gpu(s)$(gputime)\n")
end  

############################# Only Maximum Case ###################################
for i in 1:10 
    export_matrix = false
    global cputime = @elapsed global lodcmax = scan(Y, G, maf_threshold=0.05, export_matrix=export_matrix, usegpu=false, lod_or_pval="lod", timing_file=timing_file)
    global gputime = CUDA.@elapsed global lodgmax = scan(Y, G,  maf_threshold=0.05, usegpu=true, export_matrix=export_matrix, timing_file=timing_file);
end

@debug begin 
    "Results of compraring CPU full lod score to GPU full lod score matrix: $(sum(isapprox.(lodcmax, lodgmax, atol=0.001)) == prod(size(lodcmax)))"
end

open(timing_file, "a") do io
    write(io, "#$(now()),n:$n m:$m p:$p,LiteQTL Only Max,nthreads: $(nthreads()),cpu(s):$(cputime),gpu(s)$(gputime)\n")
end   

# liteqtllod = scan(Y, G, maf_threshold=0.00, export_matrix=true, usegpu=false, lod_or_pval="lod", timing_file=timing_file)

# writedlm("julia-scan-lod-result-cpu.csv", liteqtllod, ',')