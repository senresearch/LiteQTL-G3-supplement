

using CSV
using DataFrames

function check_results(rqtl_result::DataFrame, julia_result::DataFrame)
    #remove row name column. 
    rqtl_float = rqtl_result |> Matrix |> transpose |> collect
    rqtl_max = find_max_idx_value(rqtl_float)

    julia_max = Matrix(julia_result)
    for i in size(julia_max)[1]
        # compare value first 
        if !isapprox(rqtl_max[i,2], julia_max[i,2], atol=1e-5)
            error("Scan result does not agree at [$i,$j], rqtl: $(rqtl_max[i,2]), julia: $(julia_max[i,2])")
            return;
        end

        # then compare index 
        if rqtl_max[i, 1] != julia_max[i,1]
            print("Scan result index is different, rqtl result is at index: $(rqtl_max[i, 1]), julia result is at index: $(julia_max[i,1])")
        end

    end


    return "Scan result agrees. "
end


# function check_results(rqtl_result::AbstractString, julia_result::AbstractString)
#     #remove row name column. 
#     rqtl_float = rqtl_result |> Matrix |> transpose
#     julia_float = Matrix(julia_result)
#     for j in size(rqtl_float)[2]
#         for i in size(rqtl_float)[1]
#             if !isapprox(rqtl_float[i,j], julia_float[i,j], atol=1e-5)
#                 error("Scan result does not agree at [$i,$j], rqtl: $(rqtl_float[i,j]), julia: $(julia_float[i,j])")
#                 return;
#             end
#         end
#     end
#     return "Scan result agrees. "
# end
function find_max_idx_value(lod::AbstractArray{<:Real,2})
    max_array = Array{typeof(lod[1,1]),2}(undef, size(lod)[1], 2)
    Threads.@threads for i in 1:size(lod)[1]
        temp = lod[i, 1]
        idx = 1
        for j in 2:size(lod)[2]
            if temp < lod[i,j]
                temp = lod[i,j]
                idx = j
            end
        end
        max_array[i,1] = idx
        max_array[i,2] = temp
    end
    return max_array
end


for datatype in [Float64, Float32]
    for dataset in ["spleen", "hippo"]
        rqtl_result_file = joinpath(Base.@__DIR__, "..", "data", "results", dataset*"_rqtl_lod_score.csv")
        julia_result_file = joinpath(Base.@__DIR__, "..", "data", "results", string(datatype) * dataset*"_lmgpu_output.csv")
        
        rqtl_result =  CSV.read(rqtl_result_file, datarow=2)
        julia_result = CSV.read(julia_result_file, datarow=1)

        check_results(rqtl_result, julia_result)
    end
end


