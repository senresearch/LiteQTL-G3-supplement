all: installlib installcuda getdata cleandata rscan juliascan gpuscan compare gemm
nogpu: installlib getdata cleandata rscan juliascan compare

.PHONY: all



installlib:
	Rscript rlibdep.R 2>&1 | tee ../log/installlib.log
	julia juliadep.jl 2>&1 | tee -a ../log/installlib.log

installcuda:
	julia installcuda.jl 2>&1 | tee ../log/installcuda.log

getdata:
	julia download-data.jl 2>&1 | tee ../log/download-data.log

cleandata:
	Rscript clean-data.R 2>&1 | tee ../log/clean-data.log

rscan:
	# running the scan will also generate geno prob and pheno no missing data file. 
	Rscript rqtl2scan.R 2>&1 | tee ../log/rscan.log

juliascan:
	JULIA_NUM_THREADS=16 julia julia_genome_scan.jl 2>&1 | tee ../log/juliascan.log
gpuscan: 
	JULIA_NUM_THREADS=16 julia julia_genome_scan_gpu.jl 2>&1 | tee ../log/gpuscan.log
compare:
	julia check_results.jl 2>&1 | tee ../log/compare.log
gemm: 
	JULIA_NUM_THREADS=16 julia dgemm-test.jl 2>&1 | tee ../log/gemm.log
	julia generate_figure.jl 2>&1 | tee -a ../log/gemm.log

 
