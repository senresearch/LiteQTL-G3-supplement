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
    writedlm(joinpath(Base.@__DIR__, "..", "data", "results", string(datatype) * output_file), lod, ',')

end

spl_geno_file = joinpath(@__DIR__, "..", "data", "processed", "bxd-genoprob_hippo.csv")
spl_pheno_file = joinpath(@__DIR__, "..", "data","processed", "hippo-pheno-nomissing.csv")

spl_output_file = "lmgpu_hippo_output.csv"
spl_export_matrix = false

main_scan(spl_geno_file, spl_pheno_file, spl_output_file, spl_export_matrix, Float64)
main_scan(spl_geno_file, spl_pheno_file, spl_output_file, spl_export_matrix, Float32)


hip_geno_file = joinpath(@__DIR__, "..", "data", "processed", "bxd-genoprob_hippo.csv")
hip_pheno_file = joinpath(@__DIR__, "..", "data","processed", "hippo-pheno-nomissing.csv")

hip_output_file = "lmgpu_hippo_output.csv"
# don't set hip_export_matrix to true because the result matrix can take too much space on disk. 
hip_export_matrix = false

main_scan(hip_geno_file, hip_pheno_file, hip_output_file, hip_export_matrix, Float64)
main_scan(hip_geno_file, hip_pheno_file, hip_output_file, hip_export_matrix, Float32)
