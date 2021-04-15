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

open(timing_file, "a") do io
    write(io, "Time,Device,data_transfer_time,compute_time,result_reorg_time,pval_time,elapsed_total\n")
end 


############################# Full Matrix Case ###################################
cputime = 0 
gputime = 0
for i in 1:10
    export_matrix = true
    global putime = @elapsed lodcfull = scan(Y, G, maf_threshold=0.0, export_matrix=export_matrix, usegpu=false, lod_or_pval="lod", timing_file=timing_file)
    global gputime = CUDA.@elapsed lodgfull = scan(Y, G,  maf_threshold=0.0, usegpu=true, export_matrix=export_matrix, timing_file=timing_file);
    # gputime = CUDA.@elapsed lodg = scan(Y, G, maf_threshold=0.0, usegpu=true, export_matrix=false);
end
open(timing_file, "a") do io
    write(io, "#$(now()),n:$n m:$m p:$p,LiteQTL Full Matrix,nthreads: $(nthreads()),cpu(s):$(cputime),gpu(s)$(gputime)\n")
end  
# writedlm("julia-scan-pval-result-cpu.csv", lodc, ',')

############################# Only Maximum Case ###################################
for i in 1:10 
    export_matrix = false
    global cputime = @elapsed lodcmax = scan(Y, G, maf_threshold=0.5, export_matrix=export_matrix, usegpu=false, lod_or_pval="lod", timing_file=timing_file)
    global gputime = CUDA.@elapsed lodgmax = scan(Y, G,  maf_threshold=0.5, usegpu=true, export_matrix=export_matrix, timing_file=timing_file);
end

 
open(timing_file, "a") do io
    write(io, "#$(now()),n:$n m:$m p:$p,LiteQTL Only Max,nthreads: $(nthreads()),cpu(s):$(cputime),gpu(s)$(gputime)\n\n")
end   

timing_df = CSV.read(timing_file, comment="#",missingstring="NA")
fullmat_comp_cpu = timing_df[1:2:20, :compute_time]
fullmat_dt_gpu = timing_df[2:2:21, :data_transfer_time]
fullmat_comp_gpu = timing_df[2:2:21, :compute_time]
fullmat_elapsed_cpu = timing_df[1:2:20, :elapsed_total]
fullmat_elapsed_gpu = timing_df[2:2:21, :elapsed_total]
 
maxonly_comp_cpu = timing_df[21:2:40, :compute_time]
maxonly_dt_gpu = timing_df[22:2:end, :data_transfer_time]
maxonly_comp_gpu = timing_df[22:2:end, :compute_time]
maxonly_elapsed_cpu = timing_df[21:2:40, :elapsed_total]
maxonly_elapsed_gpu = timing_df[22:2:end, :elapsed_total]

println("Full Matrix: CPU elapsed runtime $(median(fullmat_elapsed_cpu)) = 0 (data transfer) + $(median(fullmat_comp_cpu)) (compute time)")
println("Full Matrix: GPU elapsed runtime $(median(fullmat_elapsed_gpu)) = $(median(fullmat_dt_gpu)) (data transfer) + $(median(fullmat_comp_gpu)) (compute time)")
println("Max Only: CPU elapsed runtime $(median(maxonly_elapsed_cpu)) = 0 (data transfer) + $(median(maxonly_comp_cpu)) (compute time)")
println("Max Only: GPU elapsed runtime $(median(maxonly_elapsed_gpu)) = $(median(maxonly_dt_gpu)) (data transfer) + $(median(maxonly_comp_gpu)) (compute time)")
