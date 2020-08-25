using PyCall
pygui(:tk)
using PyPlot
using DataFrames, CSV

pwd()

dataFileName = "dgemm_timing_12gb.csv"
df = CSV.read(dataFileName, header = ["m","n","p","gb","cpu","gpu","speedup"],  datarow = 2, limit =23);

describe(df)

df



kVec = ["k", "k", "k", "k", "k", "k", "k", "k", "k", "k", "k", "k", "k", "k", "k", "k", "k", "k", "k", "k", "k"];

mDim =(collect(df[:m]./1000))
mDim = round.(mDim,digits = 0)
mDim  = convert(Vector{Int64}, mDim)
mDim = string.(mDim).*kVec

nDim =(collect(df[:n]./1000))
nDim = round.(nDim,digits = 0)
nDim  = convert(Vector{Int64}, nDim)
nDim = string.(nDim).*kVec

pDim =(collect(df[:p]./1000))
pDim = round.(pDim,digits = 0)
pDim  = convert(Vector{Int64}, pDim)
pDim = string.(pDim).*kVec


#mDim = Int(round(collect(df[:m])./1000, digits = 0))

xname = string.(mDim,", ",nDim,", ",pDim)

speedup = log2.(collect(df[:speedup]));

xname[1] = string(xname[1]*"\n test")



fig = figure("pyplot_barplot",figsize=(10,8), dpi= 300)
b = bar(xname,speedup,color="#0f87bf",align="center",alpha=0.4);
axis("tight")
title("")
grid("off")

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



fig = figure("pyplot_barplot",figsize=(10,5), dpi= 600)
b = bar(xname,speedup,color="#0f87bf",align="center",alpha=0.4);
axis("tight")
title("")
grid("off")
xlabel("Dimensions of matrix: m, n, p")
ylabel("Speedup  log[CPU time] / [GPU time]")
fig[:autofmt_xdate](bottom=0.25,rotation=30,ha="right")



#savefig("/home/faragegr/Project/smartbandactivity2018/develop/test/temp_images/fig1.png", dpi= 600)
#gcf() # Needed for IJulia to plot inline


