# Need to install matplotlib  ` pip install matplotlib`
# generate figure for dgemm 
using DataFrames, CSV

using PyCall
pygui(:tk)
using PyPlot



gemm = CSV.read("../data/gemm-timing/dgemm-result.csv", missingstring="NA", header=0, comment="#")
header = [:m, :n, :p, :req_mem_gb, :avail_mem_gb, :cpu, :gpu, :speedup]
names!(gemm, header)

gemm_no_missing = dropmissing(gemm)

# get the data point that uses more than 8 gb of memory in total ( (m*n + n*p + m*p) * sizeof(Float64) )
gemm_no_missing[gemm_no_missing.req_mem_gb .> 8, :]

# get the ones that have a speedup of over 1.7 
# gemm_no_missing[gemm_no_missing.speedup .> 1.7, :]
selected = gemm_no_missing[gemm_no_missing.req_mem_gb .> 9, :]

function num2k_string(number::Int)
    if number > 1000
        return string(number)[1:end-3]*"k"
    else 
        return string(number)
    end
end

mDim = num2k_string.(selected[:m])
nDim = num2k_string.(selected[:n])
pDim = num2k_string.(selected[:p])

xname = string..(mDim,", ",nDim,", ",pDim)
speedup = log2.(collect(selected[:speedup]));

