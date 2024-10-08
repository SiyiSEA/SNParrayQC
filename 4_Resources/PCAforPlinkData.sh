#!/bin/bash

PCAforPlinkData () {
    bfile=$1
    prefix=$2
    nSD=$3

    # remove the high LD region (optional)
    ${PLINK}/plink \
        --bfile ${bfile} \
        --indep 50 5 1.5 \
        --exclude ${HighLD} \
        --out ${RESULTSDIR}/PCAVariants/${prefix}.ld

    ${PLINK}/plink \
        --bfile ${bfile} \
        --extract ${RESULTSDIR}/PCAVariants/${prefix}.ld.prune.in \
        --make-bed \
        --out ${RESULTSDIR}/PCAVariants/${prefix}.ld.prune

    ${GCTA}/gcta-1.94.1 \
        --bfile ${RESULTSDIR}/PCAVariants/${prefix}.ld.prune \
        --make-grm-bin \
        --autosome \
        --out ${RESULTSDIR}/PCAVariants/${prefix}.imqc

    ${GCTA}/gcta-1.94.1 \
        --grm ${RESULTSDIR}/PCAVariants/${prefix}.imqc \
        --pca \
        --out ${RESULTSDIR}/PCAVariants/${prefix}.imqc.pca

    # plot PCs to identify outliers
    Rscript ${RESOURCEDIR}/plotPCs.r \
        ${RESULTSDIR}/PCAVariants/${prefix}.imqc.pca.eigenvec \
        $nSD

    nOutlierPC1=$(grep -c "\<PC1\>" ${RESULTSDIR}/PCAVariants/${prefix}_OutliersFromPC_${nSD}SDfromMean.txt)
    echo "There are" ${nOutlierPC1} "outliers for " ${prefix} " on PC1."
    awk '$4=="PC1" {print $2,$3}' ${RESULTSDIR}/PCAVariants/${prefix}_OutliersFromPC_${nSD}SDfromMean.txt > OutlierPC1.txt
    nOutlierPC2=$(grep -c "\<PC2\>" ${RESULTSDIR}/PCAVariants/${prefix}_OutliersFromPC_${nSD}SDfromMean.txt)
    echo "There are" ${nOutlierPC2} "outliers for " ${prefix} " on PC2."
    awk '$4=="PC2" {print $2,$3}' ${RESULTSDIR}/PCAVariants/${prefix}_OutliersFromPC_${nSD}SDfromMean.txt > OutlierPC2.txt
    sort OutlierPC1.txt OutlierPC2.txt | uniq > OutlierQC.txt
    rm ${RESULTSDIR}/PCAVariants/${prefix}.ld* 
    #rm ${RESULTSDIR}/PCAVariants/${prefix}.imqc*

}
