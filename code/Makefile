all: getdata installrlib cleandata rscan #compare

.PHONY: all

getdata:
	julia download-data.jl

installrlib:
	Rscript rlibdep.R

cleandata:
	Rscript clean-spleen-data.R
	Rscript clean-hippo-data.R

rscan:
	Rscript Rqtl2scan_spleen.R
	# Rscript Rqtl2scan_hippo.R

gemmascan:

# juliascan:


# compare:
