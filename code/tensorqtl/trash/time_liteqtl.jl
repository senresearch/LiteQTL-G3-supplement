using LiteQTL 
using CSV 
using DelimitedFiles
using DataFrames
using CUDA

pheno_file = "../../data/tensorqtldata/cleanpheno.csv"
geno_file = "../../data/tensorqtldata/chr9.csv"
export_matrix = false

pheno = LiteQTL.get_pheno_data(pheno_file, Float64, transposed=true)
# 10.865656 seconds (47.21 M allocations: 1.625 GiB, 12.64% gc time)
geno = CSV.read(geno_file, DataFrame)
# 109.905806 seconds (41.48 M allocations: 46.943 GiB, 0.58% gc time)
# @time genomat = convert(Matrix{Float32}, geno[:, 3:end]) |> transpose |> collect
# geno = Nothing

datatype = Float64
genomat = convert(Matrix{datatype}, geno[:, 3:end]) |> transpose |> collect
pheno = convert(Matrix{datatype}, pheno)
# gc
#genobychr = Nothing; 
#geno = Nothing;
#chr1geno = Nothing;

LiteQTL.set_blas_threads(32);

G = genomat;
Y = pheno;

genomat = Nothing
pheno = Nothing 

n = size(Y,1)
m = size(Y,2)
p = size(G,2)
println("******* Indivuduals n: $n, Traits m: $m, Markers p: $p ****************");
println("Precompiling functions.")
small_Y = Y[1:100, 1:100]
small_G = G[1:100, 1:100]

LiteQTL.scan(small_Y, small_G,n;export_matrix=export_matrix);
LiteQTL.scan(small_Y, small_G,n;usegpu=true)


println("###### Getting timing: Calclating LOD takes: ")
@time lodc = LiteQTL.scan(Y, G,n;export_matrix=export_matrix);   
writedlm("julia-scan-result-cpu.csv", lodc, ',')
# lodc = Nothing;           
#330.735274 seconds (14.96 M allocations: 83.345 GiB, 0.17% gc time) 
gtime=CUDA.@elapsed lodg = LiteQTL.scan(Y, G,n;usegpu=true)
println("###### GPU running time is $gtime")
writedlm("julia-scan-result-gpu.csv", lodg, ',')
# 11.532153 seconds (13.87 M allocations: 794.160 MiB, 3.00% gc time) (Tux03, Nvidia V100)

##################################################################################
# With Covariates
println("Running with covariates.")
covarfile = "../../data/tensorqtldata/covariates.csv"
covar = CSV.read(covarfile, DataFrame)
X = convert(Matrix{datatype}, covar[:, 2:end])|> collect
small_X = X[1:100, :]
LiteQTL.scan(small_Y, small_G, small_X,n;export_matrix=export_matrix); 
LiteQTL.scan(small_Y, small_G, small_X,n;usegpu=true); 

println("Running genome scan with covariates. Calculating LOD takes: ")
ccovartime = @elapsed lodccovar = LiteQTL.scan(Y, G,X,n;export_matrix=export_matrix); 
println("###### Running genome scan with covariates on CPU. It takes $ccovartime seconds.")
writedlm("julia-scan-covar-result-cpu.csv", lodccovar, ',')

gcovartime = CUDA.@elapsed lodgcovar = LiteQTL.scan(Y, G,X,n;usegpu=true); 
println("###### Running genome scan with covariates on GPU. It takes $gcovartime seconds.")
writedlm("julia-scan-covar-result-gpu.csv", lodgcovar, ',')
# 4.45 seconds