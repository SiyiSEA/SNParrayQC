#!/bin/bash

## DESCRIPTIONS
#This script is for identifying the outliers in individual levels

## EXECUTION
# sh 02_CheckOutliers.sh

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

echo "running the Check outliers at $datapeth"
source ${datapeth}/config
source ${RESOURCEDIR}/PCAforPlinkData.sh
touch "$logfile_02"
exec > >(tee "$logfile_02") 2>&1

cd ${PROCESSDIR}/CheckOutliers || exit
cp ${PROCESSDIR}/QCData/${FILEPREFIX}_QCd_trimmed* ${PROCESSDIR}/CheckOutliers/.

# S_threshold="0.6"
# homo_threshold="5.4"

if [ -s ${FILEPREFIX}_QCd_trimmed.bk ]
then
    rm ${FILEPREFIX}_QCd_trimmed.bk
    rm ${FILEPREFIX}_QCd_trimmed.rds
fi

Rscript ${RESOURCEDIR}/Bigsnper_identify.r \
        ${PROCESSDIR}/CheckOutliers \
        ${FILEPREFIX}_QCd_trimmed.bed \
        ${S_threshold} \
        ${homo_threshold} \
        ${FILEPREFIX}
        
mv ./*.png ./*.pdf ${RESULTSDIR}/02/.
rm ${PROCESSDIR}/CheckOutliers/${FILEPREFIX}_QCd_trimmed*

# if [ -s ${FILEPREFIX}_QCd_trimmed.keep ]
# then
#     ${PLINK}/plink  --bfile ${FILEPREFIX}_QCd_trimmed \
#                     --keep ${FILEPREFIX}_QCd_trimmed.keep \
#                     --make-bed \
#                     --out ${FILEPREFIX}_QCd_trimmed_Bigsnper
#     cp ${FILEPREFIX}_QCd_trimmed_Bigsnper* ${RESULTSDIR}/02/.
#     PCAforPlinkData ${FILEPREFIX}_QCd_trimmed_Bigsnper  ${FILEPREFIX}_QCd_trimmed_Bigsnper 3
# else
#     echo "Please detect the threshold from the Bigsnper_PCA.pdf"
# fi
