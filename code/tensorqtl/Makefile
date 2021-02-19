all: download_data

.PHONY: all

download_data:
	wget -O ../../data/tensorqtldata/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.bed https://personal.broadinstitute.org/francois/geuvadis/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.bed 
	wget -O ../../data/tensorqtldata/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.bim https://personal.broadinstitute.org/francois/geuvadis/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.bim 
	wget -O ../../data/tensorqtldata/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.fam https://personal.broadinstitute.org/francois/geuvadis/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.fam 
	wget -O ../../data/tensorqtldata/GEUVADIS.445_samples.covariates.txt https://personal.broadinstitute.org/francois/geuvadis/GEUVADIS.445_samples.covariates.txt 
	wget -O ../../data/tensorqtldata/GEUVADIS.445_samples.expression.bed.gz https://personal.broadinstitute.org/francois/geuvadis/GEUVADIS.445_samples.expression.bed.gz 

unzip_pheno:
	gunzip -c ../../data/tensorqtldata/GEUVADIS.445_samples.expression.bed.gz > ../../data/tensorqtldata/GEUVADIS.445_samples.expression.bed

plink2csv:
	# download plink software
	wget https://zzz.bwh.harvard.edu/plink/dist/plink-1.07-x86_64.zip  --no-check-certificate
	unzip plink-1.07-x86_64.zip
	./plink-1.07-x86_64/plink --noweb --bfile ../../data/tensorqtldata/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.bed --recode --tab --out genotype.csv

clean_data: 
