all: download_data unzip_pheno setupvenv extractgeno clean_data time_tensorqtl time_liteqtl

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
setupvenv:
	git clone https://github.com/broadinstitute/tensorqtl
	cd tensorqtl/
	virtualenv venv 
	source venv/bin/activate 
	pip install -r install/requirements.txt .

activatevenv: 
	cd tensorqtl/
	source venv/bin/activate 
	# pip install -r install/requirements.txt .

extractgeno: 
	python extractgeno.py
	
clean_data: 
	Rscript clean_data.R
time_tensorqtl:
	# & means send the process to background
	#sh ./stay/get_it_ios.sh &  PIDIOS=$!
	#sh ./stay/get_it_mix.sh &  PIDMIX=$!
	#wait $PIDIOS
	#wait $PIDMIX 
	./hwswinfo.sh
	python time_tensorqtl.py
	
time_liteqtl:
	./hwswinfo.sh
	JULIA_NUM_THREADS=16 julia time_liteqtl.jl