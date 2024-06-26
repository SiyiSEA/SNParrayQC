#!/bin/bash

## DESCRIPTIONS
#This script takes genotype data outputted by genome studio and performs data QC
# including: removal of duplicate samples, update sample ids, perform sex check

## EXECUTION
# sh 01_QC.sh

## SOFTWARES
# king, plink, python, gcta

## REQUIRES the following variables in config file
## INPUT raw_data bfile
## OUTPUT data_QCd PCA 

source ./config
source ${RESOURCEDIR}/PCAforPlinkData.sh
touch "$logfile_01"
exec > >(tee "$logfile_01") 2>&1
cd ${PROCESSDIR}/QCData || exit

echo "PCA the raw data---------------------------------------------------------------------------"
# PCA the raw data -- pass
PCAforPlinkData ${RAWDATADIR}/${FILEPREFIX} ${FILEPREFIX} 2

echo "Filter on Sample-level: Chech the relatedness and duplications-----------------------------------------------------"
Method -1: identify duplication or related individuals or monozygotic twins
"$KINGPATH"/king \
    -b "${RAWDATADIR}"/"${FILEPREFIX}".bed \
    --related \
    --degree 2 \
    --prefix ${PROCESSDIR}/QCData/related

if [ -s ${PROCESSDIR}/QCData/related.kin ]
then
    awk '{if ($15 == "Dup/MZTwin") print $1, $4}' ${PROCESSDIR}/QCData/related.kin > ${PROCESSDIR}/QCData/duplicationSample.tmp
    num_dup=$(wc -l ${PROCESSDIR}/QCData/duplicationSample.tmp)
    if [ $num_dup -gt 1  ]
    then
        echo "There are " $num_dup "sample(s) is/are related."
        ${PLINK}/plink --bfile ${RAWDATADIR}/${FILEPREFIX} --missing --out ${PROCESSDIR}/QCData/MissingDuplicateSamples
        python ${SCRIPTDIR}/utilitys/ExcludeDuplicates.py${PROCESSDIR}/QCData/related.kin MissingDuplicateSamples.imiss ${PROCESSDIR}/QCData/dupsToExclude.txt
        ${PLINK}/plink --bfile ${RAWDATADIR}/${FILEPREFIX} --remove ${PROCESSDIR}/QCData/dupsToExclude.txt --make-bed --out ${PROCESSDIR}/QCData/${FILEPREFIX}_update_1
    fi
else
    echo "No duplication or related individuals."
    cp ${RAWDATADIR}/${FILEPREFIX}.bed ${FILEPREFIX}_update_1.bed
    cp ${RAWDATADIR}/${FILEPREFIX}.bim ${FILEPREFIX}_update_1.bim
    cp ${RAWDATADIR}/${FILEPREFIX}.fam ${FILEPREFIX}_update_1.fam
fi

# Method -2: identify duplicates or MZ twins.
# "$KINGPATH"/king -b "${RAWDATADIR}"/"${FILEPREFIX}".bed --duplicate --prefix ${PROCESSDIR}/QCData/duplication
# if [ -s QCoutput/king.con ]
# then 

#     cut -f 1,2 QCoutput/king.con > QCoutput/duplicateSamples.tmp
#     cut -f 3,4 QCoutput/king.con >> QCoutput/duplicateSamples.tmp
#     sort QCoutput/duplicateSamples.tmp | uniq > QCoutput/duplicateSamples.txt
#     rm QCoutput/duplicateSamples.tmp

#     if [ -s QCoutput/duplicateSamples.txt ]
#     then
    
#     ## elevated missing data rate
#     ${PLINK}/plink --bfile ${RAWDATADIR}/${FILEPREFIX} --missing --out duplicateSamples

# 	## use python script to identify duplicated with greatest missingness
# 	python ${SCRIPTDIR}/utilitys/ExcludeDuplicates.py QCoutput/king.con duplicateSamples.imiss QCoutput/dupsToExclude.txt

# 	## remove duplicates
# 	${PLINK}/plink --bfile ${RAWDATADIR}/${FILEPREFIX} --remove QCoutput/dupsToExclude.txt --make-bed --out ${FILEPREFIX}_update_1
#     fi

# else
#     cp ${RAWDATADIR}/${FILEPREFIX}.bed ${FILEPREFIX}_update_1.bed
#     cp ${RAWDATADIR}/${FILEPREFIX}.bim ${FILEPREFIX}_update_1.bim
#     cp ${RAWDATADIR}/${FILEPREFIX}.fam ${FILEPREFIX}_update_1.fam

# fi

echo "Filter on Sample-level: heterozygosity, missing rate and gender check-----------------------------------------------------"
# inital filter on SNP
${PLINK}/plink --bfile ${RAWDATADIR}/${FILEPREFIX} --maf 0.35 --geno 0.05 --mind 0.02 --hwe 0.00001 --make-bed --out ${FILEPREFIX}_update_2
${PLINK}/plink --bfile ${FILEPREFIX}_update_2 --indep 50 5 1.5 --out ${FILEPREFIX}_update_2.ld

# calculate the heterozygosity rate - hetF out of range of mean-3SD and mean+3SD
${PLINK}/plink --bfile ${FILEPREFIX}_update_2 --extract ${FILEPREFIX}_update_2.ld.prune.in --het --out ${PROCESSDIR}/QCData/hetGenotypes

# elevated missing data rates - missingF > 0.01
${PLINK}/plink --bfile ${FILEPREFIX}_update_2 --missing --out ${PROCESSDIR}/QCData/missingGenotypes

# gender check - sex error
${PLINK}/plink --bfile ${RAWDATADIR}/${FILEPREFIX} --extract ${RESULTSDIR}/PCAVariants/${FILEPREFIX}.ld.prune.in --check-sex --out ${PROCESSDIR}/QCData/SexCheck
num_sex_error=$(cat ${PROCESSDIR}/QCData/SexCheck.sexcheck | grep PROBLEM | wc -l)
echo "There are " $num_sex_error " sample(s) needed to be removed due the sex error."

# generate a list of sample fail on above checking
Rscript ${RESOURCEDIR}/list_failqc.r
${PLINK}/plink --bfile ${FILEPREFIX}_update_2 --remove ${PROCESSDIR}/QCData/fail_mis_het_sex.txt --make-bed --out ${FILEPREFIX}_update_3


echo "Filter on SNP-level: duplication variants-----------------------------------------------------"
# method-1
${PLINK}/plink --bfile ${FILEPREFIX}_update_3 --list-duplicate-vars ids-only suppress-first --out ${FILEPREFIX}_update_3
num_dup_var=$(wc -l ${FILEPREFIX}_update_3.dupvar)
echo "There are " $num_dup_var " variant(s) are duplicated."
${PLINK}/plink --bfile ${FILEPREFIX}_update_3 --exclude ${FILEPREFIX}_update_3.dupvar --make-bed --out ${FILEPREFIX}_update_4

# # method-2
# ## Variants Duplication QC - remove variants at the same position (i.e. triallelic)
# awk '{if ($1 != 0) print $1":"$4}' ${FILEPREFIX}_update_3.bim > pos.tmp
# sort pos.tmp | uniq -d > dupLocs.txt

# awk -F ":" '{print $1,$2-1,$2,"set1", "set2"}' dupLocs.txt > positionsExclude.txt
# ${PLINK}/plink --bfile ${FILEPREFIX}_update_3 --exclude range positionsExclude.txt --make-bed --out ${FILEPREFIX}_update_4

# rm QCoutput/pos.tmp
# rm QCoutput/dupLocs.txt


echo "Filter on SNP-level: missing rate of variants-----------------------------------------------------"
# # Variants missing call rate QC
# previous QC step for calculating sample missing call rate - .lmiss
${PLINK}/plink --bfile ${RAWDATADIR}/${FILEPREFIX} --missing --out rawMissing

echo "Filter on SNP-level: MAF-----------------------------------------------------"
## Variants Minor allele frequencies
# the maf threshold could be 0.01
${PLINK}/plink --bfile ${RAWDATADIR}/${FILEPREFIX} --freq --out rawVariantFreq

## filter sample and variant missingness, HWE, rare variants and exclude variants with no position
${PLINK}/plink --bfile ${FILEPREFIX}_update_4 --maf 0.01 --hwe 0.000001 --mind 0.02 --geno 0.05 --make-bed --out ${FILEPREFIX}_QCd

echo "PCA the QCd data---------------------------------------------------------------------------"
PCAforPlinkData ${FILEPREFIX}_QCd ${FILEPREFIX}_QCd 2

# ## clean up intermediate files but keep log files
Rscript ${RESOURCEDIR}/QCreport01.rmd
rm ${FILEPREFIX}_update_*.b*
rm ${FILEPREFIX}_update_*.fam