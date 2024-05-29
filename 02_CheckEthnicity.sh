#!/bin/sh

## This script determines sample ethnicity by comparing to 1000G super populations

## EXECUTION
# sh SNPArray/preprocessing/2_CheckEthnicity.sh
# where 
# script needs to be executed from <git repo>/array/

## REQUIRES the following variables in config file
# RAWDATADIR, FILEPREFIX, KGG

## REQUIRES the following software
# plink,  gcta

## INPUT
# ${FILEPREFIX}_QCd # binary plink files following prelim QC

## OUTPUT
# merge1KG/${FILEPREFIX}_mergedw1000G # variants merged with 100 genomes and filtered to common, shared variants
# merge1KG/${FILEPREFIX}_mergedw1000G.pca # pca for sample and 1000 genome combined

module load R/4.2.1-foss-2022a

cd ${PROCESSDIR} || exit
mkdir -p merge1KG

# liftover to fit the 1000G
# convert the hg17.bim file into hg18.BED file
awk '{print "chr"$1, "\t", $4-1, "\t", $4, "\t", $2}' ${FILEPREFIX}_QCd.bim > hg17.BED
#wget https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg17/liftOver/hg17ToHg38.over.chain.gz
${LIFTOVER} hg17.BED hg17ToHg38.over.chain.gz Hg38.BED unMappedHg38


# check if you have althernative chr
awk '{print $1}' Hg38.BED | sort -u
awk '{OFS="\t"; print substr($1,4), $4, 0, $3}' Hg38.BED > Hg38.bim
awk '{OFS="\t"; print $4, $3}' Hg38.BED > Hg38build.txt
${PLINK}/plink --bfile ${FILEPREFIX}_QCd --update-map Hg38build.txt --make-bed --out Hg38_update


# change variant ids to chr:bp
awk '{if ($1 != 0) print $2,"chr"$1":"$4}' Hg38_update.bim > updateTo1KGFormat.txt
${PLINK}/plink --bfile Hg38_update --update-name updateTo1KGFormat.txt --make-bed --out ${FILEPREFIX}_QCd_1kgIDs

# first merge with 1000 genomes and filter variants to those in common
# need to test initially in case of error with triallelic variants
${PLINK}/plink --bfile ${FILEPREFIX}_QCd_1kgIDs --bmerge ${KGG}/1000G_gr38_maffilt.bed ${KGG}/1000G_gr38_maffilt.bim ${KGG}/1000G_gr38_maffilt.fam --maf 0.2 --geno 0.05 --make-bed --out merge1KG/mergedw1000G_test

## issue with variants at same position but different alleles - exclude these
${PLINK}/plink --bfile ${FILEPREFIX}_QCd_1kgIDs --exclude merge1KG/mergedw1000G_test-merge.missnp --make-bed --out merge1KG/${FILEPREFIX}_1kgIDs_forMerge

${PLINK}/plink --bfile merge1KG/${FILEPREFIX}_1kgIDs_forMerge --bmerge ${KGG}/1000G_gr38_maffilt.bed ${KGG}/1000G_gr38_maffilt.bim ${KGG}/1000G_gr38_maffilt.fam --maf 0.2 --geno 0.05 --make-bed --out merge1KG/${FILEPREFIX}_mergedw1000G


# LD prune
${PLINK}/plink --bfile merge1KG/${FILEPREFIX}_mergedw1000G --indep 50 5 1.5 --out merge1KG/${FILEPREFIX}_mergedw1000G.ld
${PLINK}/plink --bfile merge1KG/${FILEPREFIX}_mergedw1000G --extract merge1KG/${FILEPREFIX}_mergedw1000G.ld.prune.in --make-bed --out merge1KG/${FILEPREFIX}_mergedw1000G.ld.prune

rm merge1KG/${FILEPREFIX}_1kgIDs_forMerge*
rm merge1KG/mergedw1000G_test*

# use GCTA to calc PCs
${GCTA}/gcta-1.94.1 --bfile merge1KG/${FILEPREFIX}_mergedw1000G.ld.prune --make-grm-bin --autosome --out merge1KG/${FILEPREFIX}_mergedw1000G
${GCTA}/gcta-1.94.1 --grm merge1KG/${FILEPREFIX}_mergedw1000G --pca --out merge1KG/${FILEPREFIX}_mergedw1000G.pca

rm merge1KG/${FILEPREFIX}_mergedw1000G*grm*

# plot PCs
Rscript ${SCRIPTDIR}/4_Resources/plotEthnicity.r ${SCRIPTDIR}/3_Results/ ${PROCESSDIR}/merge1KG/${FILEPREFIX} ${KGG} 


