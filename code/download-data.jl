bxd_geno = "http://genenetwork.org/api/v_pre1/genotypes/BXD.geno"
spleen_pheno="http://datafiles.genenetwork.org/download/GN283/GN283_MeanDataAnnotated_rev081815.txt"
hippo_pheno="http://datafiles.genenetwork.org/download/GN206/GN206_MeanDataAnnotated_rev081815.txt"

download(bxd_geno, joinpath(@__DIR__,"..", "data", "raw", "bxd.geno"))
download(spleen_pheno, joinpath(@__DIR__,"..", "data", "raw", "bxdspleen.txt"))
download(hippo_pheno, joinpath(@__DIR__,"..", "data", "raw", "bxdhippo.txt"))
