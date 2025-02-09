#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=03CheckEthnicity.o
#SBATCH --error=03CheckEthnicity.e
#SBATCH --job-name=QC03CheckEthnicity

## This script determines sample ethnicity by comparing to 1000G super populations

## EXECUTION
# sh SNPArray/preprocessing/03_CheckEthnicity.sh
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

echo "checking the arguments for config file----------------------------------------------------------------------------"
datapeth=$1

if [ -z "$1" ]
then
        echo "No argument supplied"
        echo "Please input the paht of the data folder as the first argument"
		exit 1 # fail
fi

echo "running the Check Ethnicity at $datapeth"
source ${datapeth}/config

mv ./03CheckEthnicity.o ${JOBSDIR}/03CheckEthnicity.o
mv ./03CheckEthnicity.e ${JOBSDIR}/03CheckEthnicity.e

touch "$logfile_03b"
source ${RESOURCEDIR}/PCAforPlinkData.sh
exec > >(tee "$logfile_03b") 2>&1
cd ${PROCESSDIR}/CheckEthnicity || exit


echo "Liftover the QCd data---------------------------------------------------------------------------"
# liftover to GRCh38, since 1000G is based on the GRCh38

if [ -s ${RESULTSDIR}/02/${FILEPREFIX}_QCd_Re_trimmed.fam ]
then
    cp ${RESULTSDIR}/02/${FILEPREFIX}_QCd_Re_trimmed.fam ${PROCESSDIR}/CheckEthnicity/ToBeChecked.fam
    cp ${RESULTSDIR}/02/${FILEPREFIX}_QCd_Re_trimmed.bed ${PROCESSDIR}/CheckEthnicity/ToBeChecked.bed
    cp ${RESULTSDIR}/02/${FILEPREFIX}_QCd_Re_trimmed.bim ${PROCESSDIR}/CheckEthnicity/ToBeChecked.bim
else
    if [ -s ${RESULTSDIR}/01/${FILEPREFIX}_QCd_trimmed.fam ]
    then
        cp ${RESULTSDIR}/01/${FILEPREFIX}_QCd_trimmed.fam ${PROCESSDIR}/CheckEthnicity/ToBeChecked.fam
        cp ${RESULTSDIR}/01/${FILEPREFIX}_QCd_trimmed.bed ${PROCESSDIR}/CheckEthnicity/ToBeChecked.bed
        cp ${RESULTSDIR}/01/${FILEPREFIX}_QCd_trimmed.bim ${PROCESSDIR}/CheckEthnicity/ToBeChecked.bim
    else
        cp ${RESULTSDIR}/01/${FILEPREFIX}_QCd.fam ${PROCESSDIR}/CheckEthnicity/ToBeChecked.fam
        cp ${RESULTSDIR}/01/${FILEPREFIX}_QCd.bed ${PROCESSDIR}/CheckEthnicity/ToBeChecked.bed
        cp ${RESULTSDIR}/01/${FILEPREFIX}_QCd.bim ${PROCESSDIR}/CheckEthnicity/ToBeChecked.bim
    fi
fi



# convert the hg17.bim file into GRCh38.BED file
awk '{print "chr"$1, "\t", $4-1, "\t", $4, "\t", $2}' ToBeChecked.bim > QCd.BED
${LIFTOVER} QCd.BED "${LiftChain}" Mapped.BED unMapped 
mapped_variant=$(wc -l Mapped.BED)
total_variant=$(wc -l QCd.BED)
echo "${mapped_variant} variants have been mapped."
echo "${total_variant} variants in total."
# echo $(( $mapped_variant*100/$total_variant )) "% of variants have been liftovered successfully!"

# check if you have althernative chr
awk '{print $1}' Mapped.BED | sort -u
awk '{OFS="\t"; print $4, $3}' Mapped.BED > NewPosition.txt

${PLINK}/plink --bfile ToBeChecked \
                --update-map NewPosition.txt \
                --chr 1-23 \
                --make-bed \
                --allow-no-sex \
                --out ${RESULTSDIR}/03/${FILEPREFIX}_QCd_hg38

echo "Update the vairants ID for QCd data---------------------------------------------------------------------------"
# change variant ids to chr:bp since 1000G_gr38_maffilt is chr:bp
awk '{if ($1 != 0) print $2, "chr"$1":"$4}' ${RESULTSDIR}/03/${FILEPREFIX}_QCd_hg38.bim > updateTo1KGFormat.txt
${PLINK}/plink --bfile ${RESULTSDIR}/03/${FILEPREFIX}_QCd_hg38 \
                --update-name updateTo1KGFormat.txt \
                --allow-no-sex \
                --make-bed \
                --out ${FILEPREFIX}_QCd_1kgIDs

echo "Frist merge with 1000G in hg38---------------------------------------------------------------------------"
# first merge with 1000 genomes and filter variants to those in common
# need to test initially in case of error with triallelic variants
${PLINK}/plink --bfile ${FILEPREFIX}_QCd_1kgIDs \
                --bmerge ${Ref1000G} \
                --allow-no-sex \
                --make-bed \
                --out mergedw1000G_gr38maf

echo "Second merge with 1000G in hg38---------------------------------------------------------------------------"
## issue with variants at same position but different alleles - exclude these
${PLINK}/plink --bfile ${FILEPREFIX}_QCd_1kgIDs \
                --exclude mergedw1000G_gr38maf-merge.missnp \
                --make-bed \
                --allow-no-sex \
                --out ${FILEPREFIX}_1kgIDs_forMerge

${PLINK}/plink --bfile ${FILEPREFIX}_1kgIDs_forMerge \
                --bmerge ${KGG}/1000G_gr38_maffilt \
                --maf 0.1 \
                --geno 0.1 \
                --make-bed \
                --allow-no-sex \
                --out ${FILEPREFIX}_merged_1000G_forLD


# PCA - 2mins
echo "PCA the merged data---------------------------------------------------------------------------"
PCAforPlinkData ${FILEPREFIX}_merged_1000G_forLD ${FILEPREFIX}_merged_1000G_forLD 3

# plot PCs - take 11mins
echo "Plot the PCs---------------------------------------------------------------------------"
Rscript ${SCRIPTDIR}/4_Resources/plotEthnicity.r ${RESULTSDIR}/03 ${RESULTSDIR}/PCAVariants/${FILEPREFIX}_merged_1000G_forLD ${KGG} 

# remove redundants
rm ${FILEPREFIX}_1kgIDs_forMerge*
populations=($(cut -f3 --delim="," ${DATADIR}/3_Results/03/PredictedPopulations.csv | tail -n +2 | sort | uniq))
echo "The individuals from the data could belongs to" "${populations[@]}"
echo "If your data comes from more than one population, please check the relatedness for each population."