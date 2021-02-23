all: download_data unzip_pheno extractgeno clean_data

.PHONY: all

download_data:
	wget -O ../../data/tensorqtldata/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.bed https://personal.broadinstitute.org/francois/geuvadis/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.bed 
	wget -O ../../data/tensorqtldata/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.bim https://personal.broadinstitute.org/francois/geuvadis/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.bim 
	wget -O ../../data/tensorqtldata/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.fam https://personal.broadinstitute.org/francois/geuvadis/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.fam 
	wget -O ../../data/tensorqtldata/GEUVADIS.445_samples.covariates.txt https://personal.broadinstitute.org/francois/geuvadis/GEUVADIS.445_samples.covariates.txt 
	wget -O ../../data/tensorqtldata/GEUVADIS.445_samples.expression.bed.gz https://personal.broadinstitute.org/francois/geuvadis/GEUVADIS.445_samples.expression.bed.gz 

unzip_pheno:
	gunzip -c ../../data/tensorqtldata/GEUVADIS.445_samples.expression.bed.gz > ../../data/tensorqtldata/GEUVADIS.445_samples.expression.bed

# plink2csv:
# 	# download plink software
# 	# wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip  --no-check-certificate
# 	# unzip plink*
# 	./plink --bfile ../../data/tensorqtldata/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup --recode --out ../../data/tensorqtldata/genotype
extractgeno: 
	python extractgeno.py
	
clean_data: 
	Rscript clean_data.Rscript

