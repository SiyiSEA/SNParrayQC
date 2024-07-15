#!/bin/bash

## DESCRIPTIONS
#This script is for identifying the outliers in individual levels

## EXECUTION
# sh 03_CheckOutliers.sh

## SOFTWARES
# plink, R

## REQUIRES the following variables in config file
## INPUT QCd bfile
## OUTPUT data_QCd_trimmed PCA

source ./config
source ${RESOURCEDIR}/PCAforPlinkData.sh
touch "$logfile_03"
exec > >(tee "$logfile_03") 2>&1

cd ${PROCESSDIR}/CheckOutliers || exit

cp ${PROCESSDIR}/QCData/${FILEPREFIX}_QCd* ${PROCESSDIR}/CheckOutliers/.

S_threshold="0.6"
homo_threshold="5.4"

if [ -s ${FILEPREFIX}_QCd.bk ]
then
    rm ${FILEPREFIX}_QCd.bk
    rm ${FILEPREFIX}_QCd.rds
fi

Rscript ${RESOURCEDIR}/Bigsnper_identify.r \
        ${PROCESSDIR}/CheckOutliers \
        ${FILEPREFIX}_QCd.bed \
        ${S_threshold} \
        ${homo_threshold}

if [ -s ${FILEPREFIX}_QCd.keep ]
then
    ${PLINK}/plink  --bfile ${FILEPREFIX}_QCd \
                    --keep ${FILEPREFIX}_QCd.keep \
                    --make-bed \
                    --out ${FILEPREFIX}_QCd_trimmed
    cp ${FILEPREFIX}_QCd_trimmed.* ${RESULTSDIR}/03/.
    PCAforPlinkData ${FILEPREFIX}_QCd_trimmed ${FILEPREFIX}_QCd_trimmed 3
else
    echo "Please detect the threshold from the Bigsnper_PCA.pdf"
fi
