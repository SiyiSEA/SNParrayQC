#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=/lustre/home/sww208/QC/SNParrayQC/5_JobReports/job02PLINK.o
#SBATCH --error=/lustre/home/sww208/QC/SNParrayQC/5_JobReports/job02PLINK.e
#SBATCH --job-name=JOB02

source ./config
source ${RESOURCEDIR}/PCAforPlinkData.sh
mkdir -p ${PROCESSDIR}/CheckEthnicity/1000GCompare
cd ${PROCESSDIR}/CheckEthnicity/1000GCompare || exit


Ref1000G="/lustre/home/sww208/QC/SNParrayQC/4_Resources/1000G/1000GfromPLINKWeb/all_phase3_1000G_update_4"

echo "Frist merge with 1000G in hg19---------------------------------------------------------------------------"
# first merge with 1000 genomes and filter variants to those in common
${PLINK2} --bfile ${RESULTSDIR}/02/${FILEPREFIX}_QCd_hg19 \
                --set-all-var-ids @:#_\$1_\$2 \
                --make-bed \
                --out ${FILEPREFIX}_QCd_1kgIDs

# need to test initially in case of error with triallelic variants
${PLINK}/plink --bfile ${FILEPREFIX}_QCd_1kgIDs \
                --bmerge ${Ref1000G} \
                --make-bed \
                --out mergedw1000G

echo "Second merge with 1000G in hg38---------------------------------------------------------------------------"
## issue with variants at same position but different alleles - exclude these
${PLINK}/plink --bfile ${FILEPREFIX}_QCd_1kgIDs \
                --exclude mergedw1000G-merge.missnp \
                --make-bed \
                --out ${FILEPREFIX}_1kgIDs_forMerge

${PLINK}/plink --bfile ${FILEPREFIX}_1kgIDs_forMerge \
                --bmerge ${Ref1000G} \
                --maf 0.01 \
                --geno 0.1 \
                --make-bed \
                --out merged_1000G_fromPLINK

# # remove the high LD region (optional)
# ${PLINK}/plink --bfile ${FILEPREFIX}_merged_1000G_forLD \
#                 --exclude range ${HighLD} \
#                 --make-bed \
#                 --out ${FILEPREFIX}_mergedw1000G

# PCA - 2mins
echo "PCA the merged data---------------------------------------------------------------------------"
PCAforPlinkData merged_1000G_fromPLINK merged_1000G_fromPLINK 3

# plot PCs - take 11mins
echo "Plot the PCs---------------------------------------------------------------------------"
Rscript ${SCRIPTDIR}/4_Resources/plotEthnicity.r ${PROCESSDIR}/CheckEthnicity/1000GCompare ${RESULTSDIR}/PCAVariants/merged_1000G_fromPLINK ${KGG} 