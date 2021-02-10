# This code is adapted from code/original-scripts/spleen-analysis.jl 

using LiteQTL
using DelimitedFiles

include("benchmark.jl")

function main_scan(geno_file::AbstractString, pheno_file::AbstractString, output_file::AbstractString, maxlod::Bool, datatype::DataType)


    LiteQTL.set_blas_threads(16);
    # Read in data.
    G = LiteQTL.get_geno_data(geno_file,datatype)
    Y = LiteQTL.get_pheno_data(pheno_file,datatype, transposed=false)
    # getting geno and pheno file size.
    n = size(Y,1)
    m = size(Y,2)
    p = size(G,2)
    println("******* Individuals n: $n, Traits m: $m, Markers p: $p ****************");
    # cpu_timing = benchmark(5, cpurun, Y, G,n,export_matrix);

    # running analysis.
    lod = LiteQTL.cpurun(Y, G,n,maxlod);
    timing = benchmark(10, LiteQTL.cpurun, Y, G,n,maxlod)
    println("CPU timing: $(timing[3]) with $datatype")
    # write output to file
    writedlm(output_file, lod, ',')

end


for datatype in [Float64, Float32]
    for dataset in ["spleen", "hippo"]
        println("Julia Genome Scan for $dataset")
        geno_file = joinpath(@__DIR__, "..", "data", "processed", dataset*"-bxd-genoprob.csv")
        pheno_file = joinpath(@__DIR__, "..", "data","processed", dataset*"-pheno-nomissing.csv")

        output_file = joinpath(Base.@__DIR__, "..", "data", "results", string(datatype) * dataset*"_LiteQTL_output.csv")
        maxlod = true # same as maxlod = true

        res = main_scan(geno_file, pheno_file, output_file, maxlod, datatype)
    end
end

