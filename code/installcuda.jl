using Pkg

packages=["CUDA"]

installed=Pkg.installed()

for p in packages 
    if haskey(installed, p)
        println("Package $p is installed. ")
    else 
        Pkg.add(p)
    end

end