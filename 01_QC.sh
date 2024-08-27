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

echo "PCA the raw data---------------------------------------------------------------------------------------------------"
PCAforPlinkData ${RAWDATADIR}/${FILEPREFIX} ${FILEPREFIX} 2

echo "Filter on Sample-level: Check the relatedness and duplications-----------------------------------------------------"
# Method -1: identify duplication or related individuals or monozygotic twins -- apply
"$KINGPATH"/king \
    -b "${RAWDATADIR}"/"${FILEPREFIX}".bed \
    --related \
    --degree 2 \
    --prefix related

if [ -s related.kin0 ]
then
    awk 'NR>1 {print $1, $2}' related.kin0 > related.tmp
    num_dup=$(wc -l related.tmp)
    if [ -s related.tmp ]
    then
        echo "There are " $num_dup "sample(s) is/are related."
        ${PLINK}/plink --bfile ${RAWDATADIR}/${FILEPREFIX} \
                        --remove related.tmp \
                        --make-bed \
                        --out ${FILEPREFIX}_update_1
    fi
else
    echo "No related individuals."
    cp ${RAWDATADIR}/${FILEPREFIX}.bed ${FILEPREFIX}_update_1.bed
    cp ${RAWDATADIR}/${FILEPREFIX}.bim ${FILEPREFIX}_update_1.bim
    cp ${RAWDATADIR}/${FILEPREFIX}.fam ${FILEPREFIX}_update_1.fam
fi


# # Method -2: identify duplicates or MZ twins.
# "$KINGPATH"/king \
#     -b ${FILEPREFIX}_update_1.bed \
#     --duplicate \
#     --prefix duplication

# if [ -s duplication.con ]
# then 
#     tail -n +2 duplication.con > duplication.tmp
#     cut -f 1,2 duplication.tmp >> duplicateSamples.tmp
#     cut -f 3,4 duplication.tmp >> duplicateSamples.tmp
#     sort duplicateSamples.tmp | uniq > duplicateSamples.txt
#     rm duplicateSamples.tmp

#     if [ -s duplicateSamples.txt ]
#     then
    
#     ## elevated missing data rate
#     ${PLINK}/plink --bfile ${FILEPREFIX}_update_1 --missing --out duplicateSamples

# 	## use python script to identify duplicated with greatest missingness
# 	python ${RESOURCEDIR}/ExcludeDuplicates.py duplication.con duplicateSamples.imiss dupsToExclude.txt

# 	## remove duplicates
# 	${PLINK}/plink \
#         --bfile ${FILEPREFIX}_update_1 \
#         --remove dupsToExclude.txt \
#         --make-bed \
#         --out ${FILEPREFIX}_update_2
#     fi

# else
#     cp ${FILEPREFIX}_update_1.bed ${FILEPREFIX}_update_2.bed
#     cp ${FILEPREFIX}_update_1.bim ${FILEPREFIX}_update_2.bim
#     cp ${FILEPREFIX}_update_1.fam ${FILEPREFIX}_update_2.fam

# fi
# PCAforPlinkData ${FILEPREFIX}_update_1 ${FILEPREFIX}_update_1 3
echo "Filter on SNP-level: duplication variants--------------------------------------------------------------------------"
# method-1 -- apply
${PLINK}/plink --bfile ${FILEPREFIX}_update_1 --list-duplicate-vars ids-only suppress-first --out ${FILEPREFIX}_update_2
num_dup_var=$(wc -l ${FILEPREFIX}_update_2.dupvar)
echo "There are " $num_dup_var " variant(s) are duplicated."
${PLINK}/plink --bfile ${FILEPREFIX}_update_1 --exclude ${FILEPREFIX}_update_2.dupvar --make-bed --out ${FILEPREFIX}_update_3

# method-2
# Variants Duplication QC - remove variants at the same position (i.e. triallelic)
# awk '{if ($1 != 0) print $1":"$4}' ${FILEPREFIX}_update_2.bim > pos.tmp
# sort pos.tmp | uniq -d > dupLocs.txt
# num_dupLocs=$(wc -l dupLoc.txt)
# echo "There are " "$num_dupLocs" duplicated locus in the data.

# awk -F ":" '{print $1,$2-1,$2,"set1", "set2"}' dupLocs.txt > positionsExclude.txt
# ${PLINK}/plink --bfile ${FILEPREFIX}_update_2 \
#             --exclude range positionsExclude.txt \
#             --make-bed \
#             --out ${FILEPREFIX}_update_3


echo "Filter on SNP-level: missing rate of variants----------------------------------------------------------------------"
# Variants missing call rate QC
${PLINK}/plink --bfile ${RAWDATADIR}/${FILEPREFIX} --missing --out rawMissing

echo "Filter on SNP-level: MAF-------------------------------------------------------------------------------------------"
# Variants Minor allele frequencies
# the maf threshold could be 0.01
${PLINK}/plink --bfile ${RAWDATADIR}/${FILEPREFIX} --freq --out rawVariantFreq

# filter sample and variant missingness, HWE, rare variants and exclude variants with no position
awk '{if ($1 == 0) print $2}' ${FILEPREFIX}_update_3.bim > ${FILEPREFIX}_noLocPos.tmp
${PLINK}/plink --bfile ${FILEPREFIX}_update_3 \
                --maf 0.01 \
                --hwe 0.000001 \
                --mind 0.02 \
                --geno 0.05 \
                --exclude ${FILEPREFIX}_noLocPos.tmp \
                --make-bed \
                --out ${FILEPREFIX}_update_4

##remove variants with 3+ alleles
${PLINK}/plink --bfile ${FILEPREFIX}_update_4 \
                --biallelic-only strict \
                --make-bed \
                --out ${FILEPREFIX}_update_5



echo "Filter on Sample-level: heterozygosity, missing rate and gender check----------------------------------------------"
# LD prune
${PLINK}/plink  --bfile ${FILEPREFIX}_update_5 \
                --indep 50 5 1.5 \
                --out ${FILEPREFIX}_update_5.ld

# calculate the heterozygosity rate - hetF out of range of mean-3SD and mean+3SD
${PLINK}/plink  --bfile ${FILEPREFIX}_update_5 \
                --extract ${FILEPREFIX}_update_5.ld.prune.in \
                --het \
                --autosome \
                --out hetGenotypes

# elevated missing data rates - missingF > 0.01
${PLINK}/plink  --bfile ${FILEPREFIX}_update_5 \
                --missing \
                --out missingGenotypes

# gender check - sex error
# check may fail due to no polymorphic X chromosome loci in the data
${PLINK}/plink  --bfile ${RAWDATADIR}/${FILEPREFIX} \
                --extract ${FILEPREFIX}_update_5.ld.prune.in \
                --check-sex \
                --out SexCheck
if [ -s SexCheck.sexcheck ]
then
    num_sex_error=$(cat "SexCheck.sexcheck" | grep PROBLEM | wc -l)
    echo "There are " "$num_sex_error" " sample(s) needed to be removed due the sex error."     
else
    echo "Due to no polymorphic X chromosome loci, skip the sex check."
fi


# generate a list of sample fail on above checking and remove the failed samples
Rscript ${RESOURCEDIR}/list_failqc.r ${PROCESSDIR}/QCData
${PLINK}/plink  --bfile ${FILEPREFIX}_update_5 \
                --remove fail_mis_het_sex.txt \
                --make-bed \
                --out ${FILEPREFIX}_QCd

cp ${FILEPREFIX}_QCd.* ${RESULTSDIR}/01/.

echo "PCA the QCd data---------------------------------------------------------------------------------------------------"
PCAforPlinkData ${FILEPREFIX}_QCd ${FILEPREFIX}_QCd 3

echo "Trim the Ouliters based on the 3 sd--------------------------------------------------------------------------------"
if [ -s OutlierQC.txt ]
then
    ${PLINK}/plink  --bfile ${FILEPREFIX}_QCd \
                    --remove OutlierQC.txt \
                    --make-bed \
                    --out ${FILEPREFIX}_QCd_trimmed
    cp ${FILEPREFIX}_QCd_trimmed.* ${RESULTSDIR}/01/.
    PCAforPlinkData ${FILEPREFIX}_QCd_trimmed ${FILEPREFIX}_QCd_trimmed 3
else
    echo "There is no outliers can be identified by the 3SD"
fi
# ## clean up intermediate files but keep log files
# Rscript ${RESOURCEDIR}/QCreport01.r
rm ${FILEPREFIX}_update_*.*
