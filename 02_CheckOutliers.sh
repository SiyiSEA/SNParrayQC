#!/bin/bash

## DESCRIPTIONS
# This script is for identifying the and re-rmove the outliers at individual levels.
# This script is not necessary, only needed when the PCA plot from script 01 looks less nicer.
# Identifying outliers is based on the Bigsnper and PCA (script-1 done)
# Leave the S_threshold and homo_threshold as NA ----> Only plot the histgram and PCA plots for identify the threshold and outliers.
# There are two methods for removing the outliers
## Method - 1 Bigsnper
# Define the either S and homo parameter ----> histgram and PCA plots will be generated before and after remove the outliers based on the threshold.
# *.keep file will be generated for plink remove
## Method - 2 PCA 3SD (perfer)
#  Since the PCA function will generate the "OutlierQC.txt", remove the sample based that.

## EXECUTION
# sh 02_CheckOutliers.sh <data_path>

## SOFTWARES
# plink, R

## REQUIRES the following variables in config file
## INPUT QCd bfile
## OUTPUT data_QCd_trimmed PCA

echo "checking the arguments for config file----------------------------------------------------------------------------"
datapeth=$1

if [ -z "$1" ]
then
        echo "No argument supplied"
        echo "Please input the paht of the data folder as the first argument"
		exit 1 # fail
fi

echo "Checking outliers for $datapeth"
source ${datapeth}/config
source ${RESOURCEDIR}/PCAforPlinkData.sh
touch "$logfile_02"
exec > >(tee "$logfile_02") 2>&1

cd ${PROCESSDIR}/CheckOutliers || exit
if [ -s ${RESULTSDIR}/02/${FILEPREFIX}_QCd_Re_trimmed.fam ]
then
    cp ${RESULTSDIR}/02/${FILEPREFIX}_QCd_Re_trimmed.fam ${PROCESSDIR}/CheckOutliers/ToBeChecked.fam
    cp ${RESULTSDIR}/02/${FILEPREFIX}_QCd_Re_trimmed.bed ${PROCESSDIR}/CheckOutliers/ToBeChecked.bed
    cp ${RESULTSDIR}/02/${FILEPREFIX}_QCd_Re_trimmed.bim ${PROCESSDIR}/CheckOutliers/ToBeChecked.bim
else
    if [ -s ${RESULTSDIR}/01/${FILEPREFIX}_QCd_trimmed.fam ]
    then
        cp ${RESULTSDIR}/01/${FILEPREFIX}_QCd_trimmed.fam ${PROCESSDIR}/CheckOutliers/ToBeChecked.fam
        cp ${RESULTSDIR}/01/${FILEPREFIX}_QCd_trimmed.bed ${PROCESSDIR}/CheckOutliers/ToBeChecked.bed
        cp ${RESULTSDIR}/01/${FILEPREFIX}_QCd_trimmed.bim ${PROCESSDIR}/CheckOutliers/ToBeChecked.bim
    else
        cp ${RESULTSDIR}/01/${FILEPREFIX}_QCd.fam ${PROCESSDIR}/CheckOutliers/ToBeChecked.fam
        cp ${RESULTSDIR}/01/${FILEPREFIX}_QCd.bed ${PROCESSDIR}/CheckOutliers/ToBeChecked.bed
        cp ${RESULTSDIR}/01/${FILEPREFIX}_QCd.bim ${PROCESSDIR}/CheckOutliers/ToBeChecked.bim
    fi
fi

# S_threshold="0.6"
# homo_threshold="5.3"

if [ -s ToBeChecked.bk ]
then
    rm ToBeChecked.bk
    rm ToBeChecked.rds
fi

echo "identifying the outliers based on the Bigsnper----------------------------------------------------------------"
Rscript ${RESOURCEDIR}/Bigsnper_identify.r \
        ${PROCESSDIR}/CheckOutliers \
        ToBeChecked.bed \
        ${FILEPREFIX} \
        ${S_threshold} \
        ${homo_threshold} \
        

mv ./*.png ./*.pdf ${RESULTSDIR}/02/.

# Method - 1 for removing the outliers identified by the Bigsnper package
if [ -s ToBeChecked.keep ]
then
    echo "removing the outliers based on the Bigsnper----------------------------------------------------------------"
    ${PLINK}/plink  --bfile ToBeChecked \
                    --keep ToBeChecked.keep \
                    --make-bed \
                    --out ${FILEPREFIX}_QCd_Re_trimmed

    cp ${FILEPREFIX}_QCd_Re_trimmed.bim ${RESULTSDIR}/02/.
    cp ${FILEPREFIX}_QCd_Re_trimmed.bed ${RESULTSDIR}/02/.
    cp ${FILEPREFIX}_QCd_Re_trimmed.fam ${RESULTSDIR}/02/.

    PCAforPlinkData ${FILEPREFIX}_QCd_Re_trimmed  ${FILEPREFIX}_QCd_Re_trimmed 3
else
    echo "Please detect the threshold from the Bigsnper_PCA.pdf is you select the Bigsnper to remove the outliers"
fi


# Method - 2 for removing the outliers identified by the 3SD of PCA plot
if [ -s ${PROCESSDIR}/QCData/OutlierQC.txt  ]
then
    echo "removing the outliers based on the 3SD PCA-----------------------------------------------------------------"
    ${PLINK}/plink  --bfile ToBeChecked \
                    --remove ${PROCESSDIR}/QCData/OutlierQC.txt \
                    --make-bed \
                    --allow-no-sex \
                    --out ${FILEPREFIX}_QCd_Re_trimmed

    cp ${FILEPREFIX}_QCd_Re_trimmed.* ${RESULTSDIR}/02/.
    PCAforPlinkData ${FILEPREFIX}_QCd_Re_trimmed ${FILEPREFIX}_QCd_Re_trimmed 3
else
    echo "There is no outliers can be identified by the 3SD."
fi

# rm ToBeChecked.*