using LiteQTL 
using CSV 
using DelimitedFiles
using DataFrames

pheno_file = "../../data/tensorqtldata/cleanpheno.csv"
geno_file = "../../data/tensorqtldata/chr9.csv"
export_matrix = false

@time pheno = LiteQTL.get_pheno_data(pheno_file, Float64, transposed=true)
# 10.865656 seconds (47.21 M allocations: 1.625 GiB, 12.64% gc time)
@time geno = CSV.read(geno_file, DataFrame)
# 109.905806 seconds (41.48 M allocations: 46.943 GiB, 0.58% gc time)
# @time genomat = convert(Matrix{Float32}, geno[:, 3:end]) |> transpose |> collect
# geno = Nothing

# pick out chr1 genos 
#genobychr = groupby(geno, :Column2)
#onechr = genobychr[(Column2="9",)]
# Chr1: 5618644x447; 
# Chr2: 1738829; 
# Chr3: 919842; 
# Chr4: 944044; 
# Chr5: 820554; 
# Chr6: 849929
# Chr7: 762694
# Chr8: 719376
# Chr9: 558963
# ChrX: 436393

# It takes about 60 seconds to do this step. 
datatype = Float32
@time genomat = convert(Matrix{datatype}, geno[:, 3:end]) |> transpose |> collect
pheno = convert(Matrix{datatype}, pheno)
# gc
#genobychr = Nothing; 
#geno = Nothing;
#chr1geno = Nothing;

LiteQTL.set_blas_threads(16);

G = genomat;
Y = pheno;

#genomat = Nothing
#pheno = Nothing 

n = size(Y,1)
m = size(Y,2)
p = size(G,2)
println("******* Indivuduals n: $n, Traits m: $m, Markers p: $p ****************");
println("Precompiling functions.")
lodc = LiteQTL.scan(Y, G,n;export_matrix=export_matrix);
lodg = LiteQTL.scan(Y, G,n;usegpu=true)
lodc = Nothing; lodg = Nothing;

println("Calclating LOD takes: ")
@time lodc = LiteQTL.scan(Y, G,n;export_matrix=export_matrix);             
#330.735274 seconds (14.96 M allocations: 83.345 GiB, 0.17% gc time) 
CUDA.@elapsed lodg = LiteQTL.scan(Y, G,n;usegpu=true)
# 11.532153 seconds (13.87 M allocations: 794.160 MiB, 3.00% gc time) (Tux03, Nvidia V100)

##################################################################################
# With Covariates
println("Running with covariates.")
covarfile = "../../data/tensorqtldata/covariates.csv"
covar = CSV.read(covarfile, DataFrame)
X = convert(Matrix{datatype}, covar[:, 2:end])|> collect
covar = Nothing

println("Precompiling functions:")
lodc = LiteQTL.scan(Y, G,X,n;export_matrix=export_matrix); 
lodg = LiteQTL.scan(Y, G,X,n;usegpu=true); 
lodc = Nothing; lodg = Nothing;

println("Running genome scan with covariates. Calculating LOD takes: ")
@time lodc = LiteQTL.scan(Y, G,X,n;export_matrix=export_matrix); 
CUDA.@elapsed lodg = LiteQTL.scan(Y, G,X,n;usegpu=true); 
# 4.45 seconds