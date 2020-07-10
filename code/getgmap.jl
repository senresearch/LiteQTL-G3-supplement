using CSV 
geno_file = "../data/raw/bxd.geno"

# data starts at row 21, previous rows are comments. 
geno = readdlm(geno_file, '\t', skipstart=1)[21:end, :]
gmap = geno[:, 1:4]
writedlm("../data/processed/gmap.csv", gmap, ',')