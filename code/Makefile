all: installlib installcuda getdata cleandata rscan juliascan compare 
nogpu: installlib getdata cleandata rscan juliascan gpuscan compare gemm

.PHONY: all



installlib:
	Rscript rlibdep.R |& tee ../log/$(hostname).txt
	julia juliadep.jl |& tee -a ../log/$(hostname).txt

installcuda:
	julia installcuda.jl |& tee -a ../log/$(hostname).txt

getdata:
	julia download-data.jl |& tee -a ../log/$(hostname).txt

cleandata:
	Rscript clean-data.R |& tee -a ../log/$(hostname).txt

rscan:
	# running the scan will also generate geno prob and pheno no missing data file. 
	Rscript rqtl2scan.R |& tee -a ../log/$(hostname).txt

juliascan:
	JULIA_NUM_THREADS=16 julia julia_genome_scan.jl |& tee -a ../log/$(hostname).txt
	# runs abount 2 seconds  for max lod and float64 on my home imac 
gpuscan: 
	JULIA_NUM_THREADS=16 julia julia_genome_scan_gpu.jl |& tee -a ../log/$(hostname).txt
compare:
	julia check_results.jl |& tee -a ../log/$(hostname).txt
gemm: 
	JULIA_NUM_THREADS=16 julia dgemm-test.jl |& tee -a ../log/$(hostname).txt
	julia generate_figure.jl |& tee -a ../log/$(hostname).txt

 