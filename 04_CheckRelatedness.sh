#!/bin/bash
## because our sample is mixed ethnicity this looks funky return to after checking ethnicity.

## EXECUTION
# sh 02_CheckRelatedness.sh

## REQUIRES the following files
# ${PROCESSDIR}/QCData/${FILEPREFIX}_QCd

## REQUIRES the following software
# king, plink, R

## INPUT
# ${PROCESSDIR}/QCData/${FILEPREFIX}_QCd # binary plink files following prelim QC

## OUTPUT
# ${PROCESSDIR}/CheckRelatedness/${FILEPREFIX}_${2}_QCd_king
# ${PROCESSDIR}/CheckRelatedness/${FILEPREFIX}_${2}_QCd_ibd

source ./config
touch "$logfile_04"
exec > >(tee "$logfile_04") 2>&1
cd ${PROCESSDIR}/CheckRelatedness || exit

check_relatedness () {

   echo "Check the relatedeness for each population--------------------------------------------------------------------------"

   ${PLINK}/plink --bfile ${PROCESSDIR}/QCData/${FILEPREFIX}_QCd --keep $1 --make-bed --out ${FILEPREFIX}_${2}_QCd

   ## check for relatedness with other samples with KING
   "$KINGPATH"/king -b ${FILEPREFIX}_${2}_QCd.bed --kinship --prefix ${FILEPREFIX}_${2}_QCd_kingship

   Rscript ${SCRIPTDIR}/4_Resources/plotKinshipCoeff.r ${FILEPREFIX}_${2}_QCd_kingship.kin0 ${SCRIPTDIR}/3_Results/03

   ## check for relatedness with other samples with plink
   ${PLINK}/plink --bfile ${FILEPREFIX}_${2}_QCd --genome --mind 0.2 --out ${FILEPREFIX}_${2}_QCd_ibd

   Rscript ${SCRIPTDIR}/4_Resources/Plot_ibd.r ${FILEPREFIX}_${2}_QCd_ibd.genome ${SCRIPTDIR}/3_Results/04 $2

   rm ${FILEPREFIX}_${2}_QCd.*


}

echo "Identify the population---------------------------------------------------------------------------------------------------"
populations=($(cut -f3 --delim="," ${SCRIPTDIR}/3_Results/03/PredictedPopulations.csv | tail -n +2 | sort | uniq))

for each in ${populations[@]}
do
   grep ${each} ${SCRIPTDIR}/3_Results/03/PredictedPopulations.csv | cut -f1-2 --delim="," --output-delimiter=" " > ${each}Samples.txt
   
   
   if [[ $(wc -l <${each}Samples.txt) -ge 1 ]]
   then
      echo "Checking the relatedness on ", ${each}, "population."
      check_relatedness ${each}Samples.txt ${each}
   fi
   
done  








