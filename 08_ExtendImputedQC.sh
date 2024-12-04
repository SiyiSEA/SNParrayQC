#!/bin/bash
# this script is to assemble the imputed data
## EXECUTION
# bash ./08_ExtendImputedQC.sh

## REQUIRES the following variables in config file
# ${IMPUTEDIR}, ${RAWDATADIR}/${FILEPREFIX}, ${SCRIPTDIR}

## REQUIRES the following software
# plink1.9, plink2, GCTA, R

## INPUT
# raw data
# QCd data
# Imputed data on Michigan and Sanger
# ImputedPost data on Michigan and Sanger

## OUTPUT
# 08Extend folders in the 3_Results
# three of them are for QCdata, ImputedData and ImputedPostData
# one is PCAVariants for sorting PCA plots, outliers variants and number of variants

datapeth="/lustre/home/sww208/QC/QCDataSets/EXTEND/"
echo "running the  EXTEND IMPUTED DATA at $datapeth"
source ${datapeth}/config
mkdir -p ${RESULTSDIR}/08
logfile_08="${RESULTSDIR}/08/log08.txt"
touch "$logfile_08"
exec > >(tee "$logfile_08") 2>&1


# update the SNP ID to make sure the SNP ID of r2 file can be matched to the .bim file

echo "SNP R2" > EXTEND.all.info

cd ${DATADIR}/r2 || exit 1

for i in $(seq 1 23)
do
    echo "Combing the chr${i}.r2..."
	awk '{print $3"_"$5"_"$4, $6}' OFS='\t' chr${i}.r2 >> EXTEND.all.temp.info
done

sed 's/X/23/g' EXTEND.all.temp.info > EXTEND.all.info

wc -l EXTEND.all.info
uniq EXTEND.all.info | wc -l

# for PLINK file
${PLINK2} --bfile ${RAWDATADIR}/${FILEPREFIX} \
            --make-bed \
            --out ${FILEPREFIX}_updateID_temp \
            --set-all-var-ids @:#_\$r_\$a \
            --new-id-max-allele-len 662

${PLINK}/plink --bfile ${FILEPREFIX}_updateID_temp \
            --make-bed \
            --out ${FILEPREFIX}_updateID

# calculate the MAF
${PLINK}/plink --bfile ${FILEPREFIX}_updateID  \
                --freq \
                --out ${FILEPREFIX}_updateID_freq

rm ${FILEPREFIX}_updateID_temp*

# R script for filter the info file based on the bim file
wc -l EXTEND.all.info
wc -l ${FILEPREFIX}_updateID_freq.frq
awk 'NR==FNR { lines[$2]=1; next } $1 in lines' EXTEND.all.info EXTEND_imputed_928_sex_updateID.bim | wc -l

Rscript ${RESOURCEDIR}/EXTENDMatchR2MAF.R

# filter the plink file based on the overlapped variants
${PLINK}/plink --bfile ${FILEPREFIX}_updateID  \
                --extract EXTEND_variants_keep.txt \
                --make-bed \
                --out ${FILEPREFIX}_final

cp ${FILEPREFIX}_final.bim /lustre/home/sww208/GoDMC/DataSetGoDMC/EXTEND/input_data/EXTEND.bim
cp ${FILEPREFIX}_final.bed /lustre/home/sww208/GoDMC/DataSetGoDMC/EXTEND/input_data/EXTEND.bed
cp ${FILEPREFIX}_final.fam /lustre/home/sww208/GoDMC/DataSetGoDMC/EXTEND/input_data/EXTEND.fam
cp EXTEND.info /lustre/home/sww208/GoDMC/DataSetGoDMC/EXTEND/input_data/EXTEND.info