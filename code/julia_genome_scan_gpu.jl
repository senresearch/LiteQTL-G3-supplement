# This code is adapted from code/original-scripts/spleen-analysis.jl 

using LiteQTL
using DelimitedFiles

include("benchmark.jl")

function main_scan(geno_file::AbstractString, pheno_file::AbstractString, output_file::AbstractString, export_matrix::Bool, datatype::DataType)


    LiteQTL.set_blas_threads(16);
    # Read in data.
    G = LiteQTL.get_geno_data(geno_file,datatype)
    Y = LiteQTL.get_pheno_data(pheno_file,datatype, transposed=false)
    # getting geno and pheno file size.
    n = size(Y,1)
    m = size(Y,2)
    p = size(G,2)
    println("******* Individuals n: $n, Traits m: $m, Markers p: $p ****************");

    # running analysis.
    lodcpu = LiteQTL.scan(Y, G;export_matrix=false,maf_threshold=0.00);
    lodgpu = LiteQTL.scan(Y, G;usegpu=true,export_matrix=false,maf_threshold=0.00);
    

    if all(isapprox.(lodgpu[:,2], lodcpu[:, 2]; atol = 1e-3))
        println("results agree")
        if all(isapprox.(lodgpu[:,1], lodcpu[:, 1]; atol = 5))
            println("index agree too. ")
        else
            
            tol = 3
            num_disagree = sum(.!isapprox.(lodgpu[:,1], lodcpu[:, 1]; atol = tol))
            println("$num_disagree number of indices don't agree. ")
            
        end
    else 
        println("error in gpu reuslt checking!")
    end

    ctiming = benchmark(10, LiteQTL.scan, Y, G,export_matrix=false,maf_threshold=0.00)
    gtiming = benchmark(10, LiteQTL.scan, Y, G,usegpu=true,export_matrix=false,maf_threshold=0.00)
    println("CPU timing: $(ctiming[3]) GPU timing: $(gtiming[3]) with $datatype, $(size(lodgpu))")
    
    # write output to file
    writedlm(output_file, lodgpu, ',')
end


for datatype in [Float64, Float32]
    for dataset in ["spleen"]#, "hippo"]
        println("Julia Genome Scan for $dataset")
        geno_file = joinpath(@__DIR__, "..", "data", "processed", dataset*"-bxd-genoprob.csv")
        pheno_file = joinpath(@__DIR__, "..", "data","processed", dataset*"-pheno-nomissing.csv")

        output_file = joinpath(Base.@__DIR__, "..", "data", "results", string(datatype) * dataset*"_LiteQTL_output.csv")
        export_matrix = false # same as maxlod = true

        res = main_scan(geno_file, pheno_file, output_file, export_matrix, datatype)
    end
end


