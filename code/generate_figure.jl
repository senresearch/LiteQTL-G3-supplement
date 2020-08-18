# Need to install matplotlib  ` pip install matplotlib`
# generate figure for dgemm 
using DataFrames, CSV

using PyCall
pygui(:tk)
using PyPlot



gemm = CSV.read("../data/gemm-timing/dgemm-result.csv", missingstring="NA", header=0, comment="#")
header = [:m, :n, :p, :req_mem_gb, :avail_mem_gb, :cpu, :gpu, :speedup]
names!(gemm, header)

gemm_no_missing = dropmissing(gemm)

# get the data point that uses more than 8 gb of memory in total ( (m*n + n*p + m*p) * sizeof(Float64) )
gemm_no_missing[gemm_no_missing.req_mem_gb .> 8, :]

# get the ones that have a speedup of over 1.7 
# gemm_no_missing[gemm_no_missing.speedup .> 1.7, :]
selected = gemm_no_missing[gemm_no_missing.req_mem_gb .> 9, :]

function num2k_string(number::Int)
    if number > 1000
        return string(number)[1:end-3]*"k"
    else 
        return string(number)
    end
end

mDim = num2k_string.(selected[:m])
nDim = num2k_string.(selected[:n])
pDim = num2k_string.(selected[:p])

xname = string.(mDim,", ",nDim,", ",pDim)
speedup = log2.(collect(selected[:speedup]));


fig = figure("pyplot_barplot",figsize=(10,8), dpi= 300)
b = bar(xname,speedup,color="#0f87bf",align="center",alpha=0.4);
axis("tight")
title("")
grid(false)

##################
#  Text Styling  #
##################
font1 = Dict("family"=>"sans serif",
    "color"=>"black",
    "weight"=>"normal",
    "size"=>16)
xlabel("Dimensions of matrix: m, n, p", fontdict=font1)
ylabel("Speedup\n  log2( [CPU time] / [GPU time] )", fontdict=font1)

#fig[:autofmt_xdate](bottom=0.25,rotation=30,ha="right")

ax2 = gca()
ax2.spines["top"].set_visible(false) # Hide the top edge of the axis
ax2.spines["right"].set_visible(false) # Hide the right edge of the axis
ax2.xaxis.set_ticks_position("bottom")
ax2.yaxis.set_ticks_position("left")
ax2.spines["left"].set_position(("axes",-0.02)) # Offset the left scale from the axis
ax2.spines["bottom"].set_position(("axes",-0.05)) # Offset the bottom scale from the axis

setp(ax2.get_xticklabels(),fontsize=14,color="black")
setp(ax2.get_yticklabels(),fontsize=14,color="black")

fig[:autofmt_xdate](bottom=0.25,rotation=50,ha="right")
#ax2.set_yticks([0., 0.5, 5.])
ax2.set_aspect(aspect = "auto")

savefig("/home/faragegr/Projects/Test/speedup.png", dpi= 300)
#gcf() # Needed for IJulia to plot inline