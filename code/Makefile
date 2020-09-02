all: installlib getdata cleandata rscan juliascan compare 
nogpu: installlib getdata cleandata rscan juliascan gpuscan compare gemm

.PHONY: all

installlib:
	Rscript rlibdep.R
	julia juliadep.jl

getdata:
	julia download-data.jl

cleandata:
	Rscript clean-data.R

rscan:
	# running the scan will also generate geno prob and pheno no missing data file. 
	Rscript rqtl2scan.R

juliascan:
	JULIA_NUM_THREADS=16 julia julia_genome_scan.jl 
	# runs abount 2 seconds  for max lod and float64 on my home imac 
gpuscan: 
	JULIA_NUM_THREADS=16 julia julia_genome_scan_gpu.jl 
compare:
	julia check_results.jl
gemm: 
	JULIA_NUM_THREADS=16 julia dgemm-test.jl
	julia generate_figure.jl

 