
# generate figure for dgemm 
using DataFrames

gemm = CSV.read("../data/gemm-timing/dgemm-result.csv", missingstring="NA", header=0, comment="#")
header = [:m, :n, :p, :req_mem_gb, :avail_mem_gb, :cpu, :gpu, :speedup]
names!(gemm, header)

gemm_no_missing = dropmissing(gemm)

# get the data point that uses more than 8 gb of memory in total ( (m*n + n*p + m*p) * sizeof(Float64) )
gemm_no_missing[gemm_no_missing.req_mem_gb .> 8, :]

# get the ones that have a speedup of over 1.7 
# gemm_no_missing[gemm_no_missing.speedup .> 1.7, :]

