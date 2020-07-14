# This code is adapted from code/original-scripts/spleen-analysis.jl 

using LMGPU
using DelimitedFiles

function main()
    # if no input args.
    geno_file = joinpath(@__DIR__, "..", "data", "processed", "bxd-genoprob_spleen.csv")
    pheno_file = joinpath(@__DIR__, "..", "data","processed", "spleen-pheno-nomissing.csv")
    export_matrix = true
    output_file = "lmgpu_spleen_output.csv"
    r_sign = false

    LMGPU.set_blas_threads(16);
    # Read in data.
    G = LMGPU.get_geno_data(geno_file)
    Y = LMGPU.get_pheno_data(pheno_file)
    # getting geno and pheno file size.
    n = size(Y,1)
    m = size(Y,2)
    p = size(G,2)
    println("******* Indivuduals n: $n, Traits m: $m, Markers p: $p ****************");
    # cpu_timing = benchmark(5, cpurun, Y, G,n,export_matrix);

    # running analysis.
    lod = LMGPU.cpurun(Y, G,n,export_matrix, r_sign);

    # write output to file
    writedlm(joinpath(Base.@__DIR__, "..", "data", "results", output_file), lod, ',')

    # TODO: generate plot?
    # return lod

end


lod = main()
