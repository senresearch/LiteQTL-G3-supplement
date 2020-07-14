using CSV

rqtl_result = CSV.read("../data/results/rqtl_lod_score_spleen.csv", datarow=2, transpose=true, threaded=true)
#remove row name column. 
rqtl_result = rqtl_result[:, 2:end]
# transform to float64
rqtl_float = convert(Array{Float64,2}, rqtl_result)

julia_result = CSV.read("../data/results/lmgpu_spleen_output.csv", datarow=1, threaded=true)
julia_float = convert(Array{Float64,2}, julia_result)

for j in size(rqtl_float)[2]
    for i in size(rqtl_float)[1]
        if !isapprox(rqtl_float[i,j], julia_float[i,j], atol=1e-5)
            @error "Scan result does not agree at [$i,$j], rqtl: $(rqtl_float[i,j]), julia: $(julia_float[i,j])"
        end
    end
end

println("Scan results agree. Done.")