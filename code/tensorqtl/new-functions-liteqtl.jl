using LiteQTL 
using CSV 
using DelimitedFiles
using DataFrames
using StatsBase

using Distributions

function filter_maf(genotype::Array{Float64, 2})
    alleles = 2
    maf_threshold = Float32(0.05)
    
    af = sum(genotype, dims=1) ./ (alleles * size(genotype, 1))
    maf = replace(x -> x > 0.5 ? 1-x : x, af)

    if maf_threshold > 0 
        mask = vec(maf .>= maf_threshold)
        if size(mask,1) != size(genotype, 2)
            error("Mask dimention does not match original matrix. Mask size: $(size(mask)), Matrix size: $(size(genotype))")
        end
        genotype = genotype[:, mask]
    end
    return genotype
end


function center_normalize(mat, dims=1)
    newmat = mat .- mean(mat, dims=dims)
    return newmat ./ sqrt.(sum(newmat.^2, dims=dims))
end

function calculate_r_tensor(geno, pheno; has_covar=false)
    # remove mean from geno and pheno
    if has_covar
        pheno = pheno.mean(dims = 1) 
        geno = geno.mean(dims = 1)
    end
    return center_normalize(pheno, 1)' * center_normalize(geno, 1)
end

function pval_calc(corr, dof)
    t = corr .* sqrt.(dof ./ (1 .- corr .^2))
    pval = 2 .* cdf(TDist(dof), .-abs.(t))
    return pval
end

function get_standardized_matrix_gpu(m::AbstractArray{<:Real,2})
    return (m .- mean(m, dims=1)) ./ std(m, corrected=false, dims=1) 
end

function get_standardized_matrix(mat::AbstractArray{<:Real,2})
    m = mat
    Threads.@threads for col in 1:size(m)[2]
        summ = 0.0f0
        rows = size(m)[1]
        for row in 1:rows
            summ += m[row, col]
        end
        mean = summ/rows
        sums = 0.0f0
        for row in 1:rows
            sums += (m[row,col] - mean)^2
        end
        std = sqrt(sums/rows)
        for row in 1:rows
            m[row,col] = (m[row,col]-mean)/std
        end
    end
    return m
end

function lod2p(lod)
    return 1 .- cdf(Chisq(1),2 .* log(10) .* lod)
end
