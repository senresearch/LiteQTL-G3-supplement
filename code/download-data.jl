include("data-loc.jl")

download(bxd_geno, joinpath(@__DIR__,"..", "data", "raw", "bxd.geno"))
download(spleen_pheno, joinpath(@__DIR__,"..", "data", "raw", "bxdspleen.txt"))
download(hippo_pheno, joinpath(@__DIR__,"..", "data", "raw", "bxdhippo.txt"))
