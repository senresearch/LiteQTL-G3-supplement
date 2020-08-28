# This code is adapted from code/original-scripts/spleen-analysis.jl 

using LMGPU
using DelimitedFiles

function main_scan(geno_file::AbstractString, pheno_file::AbstractString, output_file::AbstractString, export_matrix::Bool, datatype::DataType)


    LMGPU.set_blas_threads(16);
    # Read in data.
    G = LMGPU.get_geno_data(geno_file,datatype)
    Y = LMGPU.get_pheno_data(pheno_file,datatype, transposed=true)
    # getting geno and pheno file size.
    n = size(Y,1)
    m = size(Y,2)
    p = size(G,2)
    println("******* Indivuduals n: $n, Traits m: $m, Markers p: $p ****************");
    # cpu_timing = benchmark(5, cpurun, Y, G,n,export_matrix);

    # running analysis.
    cpu_timing = @elapsed lod = LMGPU.cpurun(Y, G,n,export_matrix);
    println("CPU timing is $(cpu_timing) with $datatype")
    # write output to file
    writedlm(output_file, lod, ',')

end

for datatype in [Float64, Float32]
    for dataset in ["spleen", "hippo"]
# for datatype in [Float64]
#     for dataset in ["hippo"]
        geno_file = joinpath(@__DIR__, "..", "data", "processed", dataset*"-bxd-genoprob.csv")
        pheno_file = joinpath(@__DIR__, "..", "data","processed", dataset*"-pheno-nomissing.csv")

        output_file = joinpath(Base.@__DIR__, "..", "data", "results", string(datatype) * dataset*"_lmgpu_output.csv")
        export_matrix = false

        main_scan(geno_file, pheno_file, output_file, export_matrix, datatype)
    end
end

