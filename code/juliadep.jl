using Pkg

lmgpu_githubloc = "https://github.com/senresearch/LMGPU.jl"

packages=["DelimitedFiles", "LMGPU"]

installed=Pkg.installed()

for p in packages 
    if haskey(installed, p)
        println("Package $p is installed. ")
    elseif p == "LMGPU"  
        Pkg.add(lmgpu_githubloc)
    else 
        Pkg.add(p)
    end

end