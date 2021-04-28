using DelimitedFiles 
using Distributions

function lod2p(lod)
    return 1 .- cdf(Chisq(1),2*log(10)*lod)
end

pyresult = CSV.read("tensorqtl-scan-pval.csv", header=true)
# jugresult = readdlm("julia-scan-covar-result-gpu.csv", ',')
# jucresult = readdlm("julia-scan-covar-result-cpu.csv", ',')

juresult = CSV.read("julia-scan-lod-result-cpu.csv", header=false)

tensorqtlpval = Matrix{Float64}(pyresult[:, 2:end])

liteqtlpval = Matrix{Float64}(juresult) |> transpose |> lod2p

println("Does results from tensorqtl and liteqtl agree: $(sum(isapprox.(tensorqtlpval, liteqtlpval, atol=0.01)) == prod(size(tensorqtlpval)))")