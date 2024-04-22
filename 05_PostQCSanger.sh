#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=/lustre/home/sww208/QC/OrganizedSNParray/PostQCSanger.o
#SBATCH --error=/lustre/home/sww208/QC/OrganizedSNParray/PostQCSanger.e
#SBATCH --job-name=PostQCSanger


## output files for use with Michegan Imputation Server
# this script is for deal with the output file from Michegan Imputation Server, So, please run the imputation on Michegan first and download the zip file for each chr into the ${IMPUTEDIR}/chrzip folder
# keep the passcode that Michegan server send to you!

## EXECUTION
# sh ./05_PostQCSanger.sh 
# or
# sbatch

## REQUIRES the following variables in config file
# ${IMPUTEDIR}, ${RAWDATADIR}/${FILEPREFIX}, ${SCRIPTDIR}

## REQUIRES the following software
# plink1.9, plink2, GCTA, R

## INPUT
# *.vcf download from the Sanger Server

## OUTPUT
# data_filtered.bim, data_filtered.fam, data_filtered.bed


## IMPORTANT


module load R/4.2.1-foss-2022a
source ./config

echo 'runing 5_PostQC.sh for Imputed data from Sanger'

cd ${IMPUTEDIR} || exit
mkdir -p ImputationOutputSanger
cd ImputationOutputSanger || exit


for i in {1..22}
do
	# If you need to convert your data, don't use any other flags besides `--out`
	if [ $i -eq 23 ]; then
		${PLINK2} --vcf X.vcf.gz --make-pgen --out data_chr${i}_dose --keep-allele-order
	else
		${PLINK2} --vcf ${i}.vcf.gz --make-pgen --out data_chr${i}_dose --keep-allele-order
	fi

	# output: data_chr2_dose.log data_chr2_dose.pgen data_chr2_dose.psam data_chr2_dose.pvar
	# .psam sample information file
	# .pvar variant information file

	# You need to make sure the prefix of `pfile` is different from `--out`
	# Rename the SNP IDs if necessary to avoid possible duplicates. `set-all-var-ids` will recode your variants to chr:pos:allele1_allele2 where # allele 1 and 2 are ascii sorted
	${PLINK2} --pfile data_chr${i}_dose --make-pgen --out data_chr${i}_filtered --set-all-var-ids chr@:#_\$1_\$2 --keep-allele-order --new-id-max-allele-len 390
	cp data_chr${i}_filtered.pvar data_chr${i}.info
	# output : data_chr2_filtered.pgen data_chr2_filtered.psam data_chr2_filtered.pvar

	# Keep SNPs with MAF>0.01 or info>0.8
	${PLINK2} --pfile data_chr${i}_filtered --extract-if-info "INFO>0.4" --maf 0.01 --make-pgen --out data_chr${i}_var_filtered --keep-allele-order

	# Remove duplicate variants
	${PLINK2} --pfile data_chr${i}_var_filtered --rm-dup exclude-mismatch --make-pgen --out data_chr${i}_filtered --keep-allele-order

	# Convert `pvar/pgen/psam` to `bim/bed/fam` format
	${PLINK2} --pfile data_chr${i}_filtered --make-bed --out data_chr${i}_filtered 
	# output: chr*_filtered.fam, .bim and bed files

	#Reformat pvar file to id maf r2 
	awk -v x='#CHROM' '$0~x,EOF''{print $3,$7}' < data_chr${i}_filtered.pvar |perl -pe 's/;/ /g' | awk '{print $1,$2,$3,$5,$6}'  |perl -pe 's/ID/ID TYPED RPAF ACAN/g' |perl -pe 's/RAF=//g'| perl -pe 's/INFO=//g' > data_chr${i}_filtered.info

	Rscript ${DATADIR}/4_Resources/infoFormatted.R \
			${IMPUTEDIR}/ImputationOutputSanger/data_chr${i}_filtered.info \
			${IMPUTEDIR}/ImputationOutputSanger/data_chr${i}_filtered_Sanger.info

done

# RefPanelAF= Allele frequency in imputation reference panel
# AN = Total number of alleles in called genotypes
# AC = Allele count in genotypes
# INFO = R^2 ?
# TYPED = Site was genotyped prior to imputation
############### vs ##############
# AF = Estimated Alternate Allele Frequency
# MAF = Estimated Minor Allele Frequency 
# R2 = Estimated Imputation Accuracy (R-square)

# clean up redundant files
rm data_chr*_dose* data_chr*_var_filtered* 

# Merge them into one dataset named data_filtered_Sanger
for i in {2..22}
do 
    echo "data_chr${i}_filtered"
done > mergefile.txt

${PLINK2} --bfile data_chr1_filtered --pmerge-list mergefile.txt --make-bed --out data_filtered_Sanger
# output: .bim, bed and fam files

# The FIDs in the fam file are set to 0.
head data_filtered_Sanger.fam
cp data_filtered_Sanger.fam data_filtered_Sanger.fam.orig
awk '{print $2,$2,$3,$4,$5,$6}' < data_filtered_Sanger.fam.orig > data_filtered_Sanger.fam
# correct the FID and IID for fam file
Rscript ${DATADIR}/4_Resources/correctFIDIID.R \
	data_filtered_Sanger.fam \
	${RAWDATADIR}/${FILEPREFIX}.fam

# Combine info files into a single file
head -n1 data_chr1_filtered_Sanger.info > data_filtered_Sanger.info

for i in {1..22}
do
    awk ' NR>1 {print $0}' < data_chr${i}_filtered_Sanger.info |cat >> data_filtered_Sanger.info
done
rm data_chr*_filtered_Sanger.info data_filtered_Sanger.fam.orig data_chr*_filtered*

# check whether the number of variants is the same
echo "The number of variants in bim file is" $(wc -l data_filtered_Sanger.bim)
echo "The number of variants in info file is" $(wc -l data_filtered_Sanger.info)


# PCA
${PLINK}/plink --bfile data_filtered_Sanger --indep 50 5 1.5 --out data_filtered_Sanger.ld
${PLINK}/plink --bfile data_filtered_Sanger --extract data_filtered_Sanger.ld.prune.in --make-bed --out data_filtered_Sanger.ld.prune

mkdir -p GCTA

${GCTA}/gcta-1.94.1 --bfile data_filtered_Sanger.ld.prune --make-grm-bin --autosome --out GCTA/data_filtered_Sanger_imqc
${GCTA}/gcta-1.94.1 --grm GCTA/data_filtered_Sanger_imqc --pca --out GCTA/data_filtered_Sanger_imqc.pca

rm data_filtered_Sanger.ld.prune.b* data_filtered_Sanger.ld.prune.fam

# plot PCs to identify outliers
Rscript /lustre/home/sww208/QC/BrainFANS/array/SNPArray/utilitys/plotPCs.r GCTA/data_filtered_Sanger_imqc.pca.eigenvec 3

