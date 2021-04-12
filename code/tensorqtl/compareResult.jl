using DelimitedFiles 
using Distributions

function lod2p(lod)
    return 1-cdf(Chisq(1),2*log(10)*lod)
end

pyresult = readdlm("tensorqtl-scan-covar.csv", ',')
# jugresult = readdlm("julia-scan-covar-result-gpu.csv", ',')
# jucresult = readdlm("julia-scan-covar-result-cpu.csv", ',')

juresult = readdlm("julia-scan-pval-result-cpu.csv", ",")

pyp = pyresult[:, 4]
# jup = lod2p.(juresult[:,2])

# julia> jup[1:10]
# 10-element Array{Float64,1}:
#  2.6511109997184867e-5
#  0.00014253024305921347
#  2.4793180719129282e-5
#  2.781509721327957e-5
#  5.521627563398468e-6
#  2.2259971643734389e-13
#  1.6991201712279747e-7
#  3.2594960064358247e-10
#  6.481299963390086e-7
#  3.567182617503217e-5

# julia> pyp[1:11]
# 11-element Array{Any,1}:
#   "pval"
#  7.452469701467925e-20
#  6.480300430135909e-8
#  4.000629114200367e-6
#  2.518314982352679e-18
#  5.400053710941933e-6
#  3.2311440125543664e-33
#  1.0041918309139447e-11
#  1.5268200471370943e-6
#  8.613801458603567e-12
#  3.399486136578877e-10