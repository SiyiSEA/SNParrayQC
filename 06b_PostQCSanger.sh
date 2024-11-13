#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=100G
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=06bPostQCSanger.o
#SBATCH --error=06bPostQCSanger.e
#SBATCH --job-name=PostQCSanger

## output files for use with Sanger Imputation Server
# this script is for deal with the output file from Sanger Imputation Server, 
# So, please download the zip file for each chr into the ${IMPUTEDIR}/

## EXECUTION
# sbatch ./06b_PostQCSanger.sh <reference panel>
# <reference panel> could be either HRC or 1000G

## REQUIRES the following variables in config file
# ${IMPUTEDIR}, ${RAWDATADIR}/${FILEPREFIX}, ${SCRIPTDIR}

## REQUIRES the following software
# plink1.9, plink2, GCTA, R

## INPUT
# *.vcf download from the Sanger Server

## OUTPUT
# data_filtered_Sanger.bim, data_filteredSanger_Sanger.fam, data_filtered_Sanger.bed, data_filtered_Sanger.info

echo "checking the arguments for config file---------------------------------------------"
datapeth=$1

if [ -z "$1" ]
then
        echo "No argument supplied"
        echo "Please input the paht of the data folder as the first argument"
		exit 1 # fail
fi

echo "running the PostQCSanger at $datapeth"
source ${datapeth}/config
touch "$logfile_06b"
exec > >(tee "$logfile_06b") 2>&1
module purge
module load R/4.2.1-foss-2022a


mv 06bPostQCSanger.o ${JOBSDIR}/06bPostQCSanger.o
mv 06bPostQCSanger.e ${JOBSDIR}/06bPostQCSanger.e


echo "checking the arguments for panel--------------------------------------------------"
cd ${IMPUTEDIR} || exit

panel=$2

if [ -z "$2" ]
then
        echo "No argument supplied"
        echo "Please input the second argument either 1000G nor HRC"
		exit 1 # fail
else 
		cd ${IMPUTEDIR}/ImputationOutputSanger${panel} || exit 1
fi

echo "runing PostQC for Imputed data from Sanger-----------------------------------------"

if [ -s X.vcf.gz ]
then
	chrNum=23
	# chr23 will fail on --make-pgen if there is no sex info
	awk '{print $1,$2,$5}' ${RAWDATADIR}/${FILEPREFIX}.fam > sex_chr.info
else
	chrNum=22
	echo "Warnning: There is no SEX CHR!!!!!!!!!!!!!!!!!!!!!!!!"
fi


for i in $(seq 1 $chrNum)
do
	# convert the vcf.gz into pgen formats
	if [[ $i -eq 23 && $panel == "1000G" ]]; 
	then
		# the imputed dataset by 1000G already contains the PAR1 and PAR2 region
		# Human chrX pseudoautosomal variant(s) appear in chrX on 1000G, needs split-par or merge-par
		${PLINK2} --vcf X.vcf.gz \
				--make-pgen \
				--update-sex sex_chr.info \
				--out data_chr${i}_dose_temp0 \
				--split-par 'hg19'
		rm sex_chr.info

	elif [[ $i -eq 23 && $panel == "HRC" ]]; 
	then
		${PLINK2} --vcf X.vcf.gz \
				--make-pgen \
				--update-sex sex.info \
				--out data_chr${i}_dose_temp0
		
	else
		${PLINK2} --vcf ${i}.vcf.gz \
				  --make-pgen \
				  --out data_chr${i}_dose_temp0 \
				  --keep-allele-order
	fi

	# some variant id is . instead of rsid, needs to update var id
	${PLINK2} --pfile data_chr${i}_dose_temp0 \
			  --make-pgen \
			  --out data_chr${i}_dose \
			  --set-all-var-ids @:#_\$1_\$2 \
			  --new-id-max-allele-len 700 # because the indel, only 1000G has

	# Keep SNPs with MAF>0.01, info>0.8 and hwe 1e-05
	# MAF ref:Performance of Genotype Imputation for Low Frequency and Rare Variants from the 1000 Genomes
	${PLINK2} --pfile data_chr${i}_dose \
			  --extract-if-info "INFO>0.8" \
			  --maf 0.01 \
			  --hwe 1e-5 \
			  --make-pgen \
		      --out data_chr${i}_filtered_temp2 \
			  --keep-allele-order


	# Remove duplicate variants
	${PLINK2} --pfile data_chr${i}_filtered_temp2 \
			--rm-dup exclude-mismatch \
			--make-pgen \
			--out data_chr${i}_filtered \
			--keep-allele-order


	# Convert `pvar/pgen/psam` to `bim/bed/fam` format
	${PLINK2} --pfile data_chr${i}_filtered \
			 --make-bed \
			 --out data_chr${i}_filtered

	# calculate the MAF
	${PLINK}/plink --bfile data_chr${i}_filtered  \
				   --freq \
				   --out data_chr${i}_filtered_freq

	# Reformat pvar file to id maf info
	# each SNP has different numbers of INFO varialbes, leading to info score of each SNP located at either 5th or 6th column
	awk -v x='#CHROM' '$0~x,EOF''{print $3,$7}' < data_chr${i}_filtered.pvar | perl -pe 's/;/ /g' | awk '{print $1,$5,$6}'  | perl -pe 's/ID/ID AC INFO/g' > data_chr${i}_filtered.info
	awk '$3!="" {print $1,$3}' OFS='\t' data_chr${i}_filtered.info > data_chr${i}_filtered_temp.info
	awk '$3==""' data_chr${i}_filtered.info >> data_chr${i}_filtered_temp.info
	# join doesn't work if two files are not in the same order, so I sort it
	sort data_chr${i}_filtered_temp.info > data_chr${i}_filtered_temp2.info
	# only need MAF in the final info file
	awk '{print $2,$5}' OFS='\t' data_chr${i}_filtered_freq.frq | sort > data_chr${i}_filtered_freq_temp.frq

	join -1 1 -2 1 data_chr${i}_filtered_freq_temp.frq data_chr${i}_filtered_temp2.info | tr -d 'INFO=' > data_chr${i}_filtered_Sanger_temp.info
	sed -e '1i ID MAF INFO' data_chr${i}_filtered_Sanger_temp.info  > data_chr${i}_filtered_Sanger_temp2.info

	sed 's/X/23/g' data_chr${i}_filtered_Sanger_temp2.info > data_chr${i}_filtered_Sanger.info

done



# Merge them into one dataset named data_filtered_Sanger
echo "Merge--------------------------------------------------------------------------------"
rm -f mergefile.txt
touch mergefile.txt
for i in $(seq 2 $chrNum)
do 
	if [ -f "data_chr${i}_filtered.bed" ]; then
		echo "data_chr${i}_filtered" >> mergefile.txt
	fi
done

${PLINK2} --bfile data_chr1_filtered \
		  --pmerge-list mergefile.txt \
		  --make-bed \
		  --out data_filtered_Sanger_temp

# correct chrX into chr23
awk '{print $2}' data_filtered_Sanger_temp.bim > oldidchrX.txt
sed 's/X/23/g' oldidchrX.txt > newidchr23.txt
paste oldidchrX.txt newidchr23.txt > updatechrid.txt

# correct the FID and IID for fam file
Rscript ${SCRIPTDIR}/4_Resources/correctFIDIID.r data_filtered_Sanger_temp.fam

${PLINK}/plink --bfile data_filtered_Sanger_temp \
			   --make-bed \
			   --update-ids updateFIDIID.txt \
			   --update-name updatechrid.txt 2 1 \
			   --out data_filtered_Sanger_temp1

# correct the sex info from the raw data
awk '{print $1,$2}' data_filtered_Sanger.fam > correctFIDIID.info
sort correctFIDIID.info > correctFIDIID.info.sort

awk '{print $1,$2,$5}' ${RAWDATADIR}/${FILEPREFIX}.fam > sex_raw.info
sort sex_raw.info > sex_raw.info.sort
join -1 1 -2 1 correctFIDIID.info.sort sex_raw.info.sort -a1 > sex.info
awk '{print $1,$2,$4}' sex.info > updatesex.info


# correct the phenotype info from the raw data
awk '{print $1,$2,$6}' ${RAWDATADIR}/${FILEPREFIX}.fam > pheno_raw.info
sort pheno_raw.info > pheno_raw.info.sort
join -1 1 -2 1 correctFIDIID.info.sort pheno_raw.info.sort -a1 > pheno.info
awk '{print $1,$2,$4}' pheno.info > updatepheno.info

${PLINK}/plink --bfile data_filtered_Sanger_temp1 \
			   --make-bed \
			   --update-sex updatesex.info \
			   --pheno updatepheno.info \
			   --out data_filtered_Sanger

# Combine info files into a single file
cp data_chr1_filtered_Sanger.info data_filtered_Sanger.info

for i in $(seq 1 $chrNum)
do
    if [ -f "data_chr${i}_filtered.bed" ]; then
		awk ' NR>1 {print $0}' < data_chr${i}_filtered_Sanger.info | cat >> data_filtered_Sanger.info
	fi
done

# clean up redundant files
rm data_filtered_Sanger_temp*
rm data_chr*temp* 
rm data_chr*_filtered_Sanger.info 
rm oldidchrX.txt newidchr23.txt
rm data_chr*_filtered.info
rm sex* pheno*

cp data_filtered_Sanger* ${RESULTSDIR}/06b/

# check whether the number of variants is the same
echo "The number of variants in bim file is" $(wc -l data_filtered_Sanger.bim)
echo "The number of variants in info file is" $(wc -l data_filtered_Sanger.info)
