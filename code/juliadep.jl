using Pkg

LiteQTL_githubloc = "https://github.com/senresearch/LiteQTL.jl"

packages=["DelimitedFiles", "Statistics", "CSV", "DataFrames", "PyCall", "PyPlot", "LiteQTL", "LinearAlgebra", "Random", "Dates"]

installed=Pkg.installed()

for p in packages 
    if haskey(installed, p)
        println("Package $p is installed. ")
    elseif p == "LiteQTL"  
        Pkg.add(url=LiteQTL_githubloc)
    else 
        Pkg.add(p)
    end`

end