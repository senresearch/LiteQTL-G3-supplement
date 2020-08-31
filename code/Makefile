all: getdata installlib cleandata rscan juliascan compare
.PHONY: all

getdata:
	julia download-data.jl

installlib:
	Rscript rlibdep.R
	julia juliadep.jl

cleandata:
	Rscript clean-data.R

rscan:
	# running the scan will also generate geno prob and pheno no missing data file. 
	Rscript rqtl2scan.R

juliascan:
	JULIA_NUM_THREADS=16 julia julia_genome_scan.jl 
	 # runs abount 2 seconds  for max lod and float64 on my home imac 
compare:
	julia check_results.jl
gemm: 
	JULIA_NUM_THREADS=16 julia dgemm-test.jl
	julia generate_figure.jl

 