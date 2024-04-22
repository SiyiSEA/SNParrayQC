#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=/lustre/home/sww208/QC/OrganizedSNParray/PostQC.o
#SBATCH --error=/lustre/home/sww208/QC/OrganizedSNParray/PostQC.e
#SBATCH --job-name=PostQC


## output files for use with Michegan Imputation Server
# this script is for deal with the output file from Michegan Imputation Server, So, please run the imputation on Michegan first and download the zip file for each chr into the ${IMPUTEDIR}/chrzip folder
# keep the passcode that Michegan server send to you!

## EXECUTION
# sh ./05_PostQC.sh 
# or
# sbatch

## REQUIRES the following variables in config file
# ${IMPUTEDIR}, ${RAWDATADIR}/${FILEPREFIX}, ${SCRIPTDIR}

## REQUIRES the following software
# plink1.9, plink2, GCTA, R

## INPUT
# chr.zip under the  ${IMPUTEDIR}/chrzip folder

## OUTPUT
# data_filtered.bim, data_filtered.fam, data_filtered.bed


## IMPORTANT


module load R/4.2.1-foss-2022a
module load p7zip
source ./config

echo 'runing 5_PostQC.sh for for imputed data from Michigan'

cd ${IMPUTEDIR} || exit
mkdir -p ImputationOutputMichigan
cd ImputationOutputMichigan || exit


# unzip the result
# function for unzip all the zip file
un7zip () {
	passcode=$1
	for chr in $(seq 1 22); do
		7z x -p"${passcode}" ./chr_$chr.zip;
	done
	7z x -p"${passcode}" ./chr_X.zip
  mv chrX.dose.vcf.gz chr23.dose.vcf.gz
}

# the passcode needs to be changed
un7zip 0h1uaFOaLEWEv2



for i in {1..23}
do
# If you need to convert your data, don't use any other flags besides `--out`
${PLINK2} --vcf ./chrzip/chr${i}.dose.vcf.gz --make-pgen --out data_chr${i}_dose --keep-allele-order

# You need to make sure the prefix of `pfile` is different from `--out`
# Rename the SNP IDs if necessary to avoid possible duplicates. `set-all-var-ids` will recode your variants to chr:pos:allele1_allele2 where # allele 1 and 2 are ascii sorted
${PLINK2} --pfile data_chr${i}_dose --make-pgen --out data_chr${i}_filtered --set-all-var-ids chr@:#_\$1_\$2 --keep-allele-order --new-id-max-allele-len 390
cp data_chr${i}_filtered.pvar data_chr${i}.info

# Keep SNPs with MAF>0.01 or info>0.8
${PLINK2} --pfile data_chr${i}_filtered --extract-if-info "R2>0.8" --maf 0.01 --make-pgen --out data_chr${i}_var_filtered --keep-allele-order

# Remove duplicate variants
${PLINK2} --pfile data_chr${i}_var_filtered --rm-dup exclude-mismatch --make-pgen --out data_chr${i}_filtered --keep-allele-order

# Convert `pvar/pgen/psam` to `bim/bed/fam` format
${PLINK2} --pfile data_chr${i}_filtered --make-bed --out data_chr${i}_filtered 

#Reformat pvar file to id maf r2 
awk -v x='#CHROM' '$0~x,EOF''{print $3,$7}' <data_chr${i}_filtered.pvar | perl -pe 's/;/ /g'|awk '{print $1,$3,$4}' |perl -pe 's/ID/ID MAF R2/g' |perl -pe 's/MAF=//g'| perl -pe 's/R2=//g' > data_chr${i}_filtered.info

done

# clean up redundant files
rm data_chr*_dose data_chr*_var_filtered* 

# Merge them into one dataset

for i in {2..23}
do 
    echo "data_chr${i}_filtered"
done > mergefile.txt

${PLINK2} --bfile data_chr1_filtered --pmerge-list mergefile.txt --make-bed --out data_filtered_Michigan

# The FIDs in the fam file are set to 0.

head data_filtered_Michigan.fam
cp data_filtered_Michigan.fam data_filtered_Michigan.fam.orig
awk '{print $2,$2,$3,$4,$5,$6}' < data_filtered_Michigan.fam.orig > data_filtered_Michigan.fam

Rscript /lustre/home/sww208/QC/OrganizedSNParray/4_Resources/correctFIDIID.R \
	data_filtered_Michigan.fam \
	${RAWDATADIR}/${FILEPREFIX}.fam


# Combine info files into a single file

head -n1 data_chr1_filtered.info > data_filtered_Michigan.info

for i in {1..23}
do
    awk ' NR>1 {print $0}' < data_chr${i}_filtered.info |cat >> data_filtered_Michigan.info
done

# check whether the number of variants is the same
echo "The number of variants in bim file is" $(wc -l data_filtered_Michigan.bim)
echo "The number of variants in info file is" $(wc -l data_filtered_Michigan.info)


# PCA
${PLINK}/plink --bfile data_filtered_Michigan --indep 50 5 1.5 --out data_filtered_Michigan_imqc

mkdir -p GCTA

${GCTA}/gcta-1.94.1 --bfile data_filtered_Michigan_imqc.ld.prune --make-grm-bin --autosome --out GCTA/data_filtered_Michigan_imqc
${GCTA}/gcta-1.94.1 --grm GCTA/data_filtered_Michigan_imqc --pca --out GCTA/data_filtered_Michigan_imqc.pca

## plot PCs to identify outliers
Rscript ${SCRIPTDIR}/SNPArray/utilitys/plotPCs.r data_filtered_Michigan_imqc.pca.eigenvec 3

# PCA
${PLINK}/plink --bfile data_filtered_Michigan --indep 50 5 1.5 --out data_filtered_Michigan.ld
${PLINK}/plink --bfile data_filtered_Michigan --extract data_filtered_Michigan.ld.prune.in --make-bed --out data_filtered_Michigan.ld.prune

mkdir -p GCTA

${GCTA}/gcta-1.94.1 --bfile data_filtered_Michigan.ld.prune --make-grm-bin --autosome --out GCTA/data_filtered_Michigan_imqc
${GCTA}/gcta-1.94.1 --grm GCTA/data_filtered_Michigan_imqc --pca --out GCTA/data_filtered_Michigan_imqc.pca

rm data_filtered_Michigan.ld.prune.b*
rm data_filtered_Michigan.ld.prune.fam

# plot PCs to identify outliers
Rscript /lustre/home/sww208/QC/BrainFANS/array/SNPArray/utilitys/plotPCs.r GCTA/data_filtered_Michigan_imqc.pca.eigenvec 3

