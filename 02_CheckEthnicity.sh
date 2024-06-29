#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=/lustre/home/sww208/QC/SNParrayQC/5_JobReports/02CheckEthnicity.o
#SBATCH --error=/lustre/home/sww208/QC/SNParrayQC/5_JobReports/02CheckEthnicity.e
#SBATCH --job-name=02CheckEthnicity

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

source ./config
touch "$logfile_02"
source ${RESOURCEDIR}/PCAforPlinkData.sh
exec > >(tee "$logfile_02") 2>&1
cd ${PROCESSDIR}/CheckEthnicity || exit

# # liftover to GRCh37/hg19, since 1000G is based on the GRCh37/hg19
# # convert the hg17.bim file into hg19.BED file
# awk '{print "chr"$1, "\t", $4-1, "\t", $4, "\t", $2}' ${RESULTSDIR}/01/${FILEPREFIX}_QCd.bim > QCd.BED
# ${LIFTOVER} QCd.BED "${LiftChain}" 1000Ghg19.BED unMapped1000Ghg19

# # check if you have althernative chr
# awk '{print $1}' 1000Ghg19.BED | sort -u
# awk '{OFS="\t"; print $4, $3}' 1000Ghg19.BED > hg19NewPosition.txt

# ${PLINK}/plink --bfile ${RESULTSDIR}/01/${FILEPREFIX}_QCd \
#                 --update-map hg19NewPosition.txt \
#                 --make-bed \
#                 --out ${RESULTSDIR}/02/${FILEPREFIX}_QCd_hg19

# # change variant ids to chr:bp
# awk '{if ($1 != 0) print $2, $1":"$4}' ${RESULTSDIR}/02/${FILEPREFIX}_QCd_hg19.bim > updateTo1KGFormat.txt
# ${PLINK}/plink --bfile ${RESULTSDIR}/02/${FILEPREFIX}_QCd_hg19 \
#                 --update-name updateTo1KGFormat.txt \
#                 --make-bed \
#                 --out ${FILEPREFIX}_QCd_1kgIDs

# update the variants ID for 1000G
awk '{if ($1 != 0) print $2, $1":"$4}' ${Ref1000G}.bim > 1000GNewPosition.txt
# how many vairants at the same position
awk '{print $2}' 1000GNewPosition.txt | sort -d | wc -l 
# needs to deal with the alternative variants before merge
${PLINK}/plink --bfile ${Ref1000G} \
                --update-name 1000GNewPosition.txt \
                --make-bed \
                --out ${Ref1000G}_update



# first merge with 1000 genomes and filter variants to those in common
# need to test initially in case of error with triallelic variants
${PLINK}/plink --bfile ${FILEPREFIX}_QCd_1kgIDs \
                --bmerge ${Ref1000G}_update \
                --maf 0.2 \
                --geno 0.05 \
                --make-bed \
                --out mergedw1000G_test

## issue with variants at same position but different alleles - exclude these
${PLINK}/plink --bfile ${FILEPREFIX}_QCd_1kgIDs \
                --exclude mergedw1000G_test-merge.missnp \
                --make-bed \
                --out ${FILEPREFIX}_1kgIDs_forMerge

# ${PLINK}/plink --bfile ${FILEPREFIX}_1kgIDs_forMerge \
#                 --bmerge ${Ref1000G} \
#                 --maf 0.2 \
#                 --geno 0.05 \
#                 --make-bed \
#                 --out ${FILEPREFIX}_merged_1000G_forLD

# # remove the high LD region
# ${PLINK}/plink --bfile ${FILEPREFIX}_merged_1000G_forLD \
#                 --exclude range ${HighLD} \
#                 --make-bed \
#                 --out ${FILEPREFIX}_mergedw1000G

# # PCA
# PCAforPlinkData ${FILEPREFIX}_mergedw1000G ${FILEPREFIX}_mergedw1000G 3

# plot PCs
# Rscript ${SCRIPTDIR}/4_Resources/plotEthnicity.r ${RESULTSDIR}/02 ${RESULTSDIR}/PCAVariants/${FILEPREFIX} ${KGG} 

# remove redundants
# rm mergedw1000G_test*