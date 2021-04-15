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
	git clone https://github.com/senresearch/tensorqtl ~/git/tensorqtl
	cd ~/git/tensorqtl/ && \
	git checkout 1760a34 && \
	virtualenv venv && \
	source venv/bin/activate && \
	pip install -r install/requirements.txt . && \
	pip install torch==1.7.1+cu101 torchvision==0.8.2+cu101 torchaudio==0.7.2 -f https://download.pytorch.org/whl/torch_stable.html 
	# Make sure the tensorqtl used is from your disk, and reflects changes you made on the source code. 
	pip install -e ~/git/tensorqtl/

activatevenv: 
	cd ~/git/tensorqtl/; \
	source venv/bin/activate

extractgeno: 
	python extractgeno.py
	
clean_data: 
	Rscript clean_data.R

compare_result_with_tensorqtl:
	python small-data-test-tensorqtl.py 
	julia small-data-test-liteqtl.jl 
	julia compareResult.jl 

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
