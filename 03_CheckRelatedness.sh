## because our sample is mixed ethnicity this looks funky return to after checking ethnicity.
## 

## EXECUTION
# sh SNPArray/preprocessing/3_CheckRelatedness.sh
# where 
# script needs to be executed from <git repo>/array/

## REQUIRES the following variables in config file
# RAWDATADIR, FILEPREFIX, 

## REQUIRES the following software
# king, plink, 

## INPUT
# ${FILEPREFIX}_QCd # binary plink files following prelim QC

## OUTPUT
# QCoutput/${FILEPREFIX}_${2}_QCd_king
# QCoutput/${FILEPREFIX}_${2}_QCd_ibd


check_relatedness () {
cd ${PROCESSDIR} || exit

${PLINK}/plink --bfile ${FILEPREFIX}_QCd --keep $1 --make-bed --out QCoutput/${FILEPREFIX}_${2}_QCd

## check for relatedness with other samples with KING
$KINGPATH/king -b QCoutput/${FILEPREFIX}_${2}_QCd.bed --kinship --prefix QCoutput/${FILEPREFIX}_${2}_QCd_king

Rscript ${SCRIPTDIR}/4_Resources/plotKinshipCoeff.r QCoutput/${FILEPREFIX}_${2}_QCd_king.kin0 ${SCRIPTDIR}/3_Results

## check for relatedness with other samples with plink
${PLINK}/plink --bfile QCoutput/${FILEPREFIX}_${2}_QCd --genome --mind 0.2 --out QCoutput/${FILEPREFIX}_${2}_QCd_ibd

#Rscript ${SCRIPTDIR}/4_Resources/Plot_ibd.R QCoutput/${FILEPREFIX}_${2}_QCd_ibd.genome ${SCRIPTDIR}/3_Results $2

rm QCoutput/${FILEPREFIX}_${2}_QCd.*


}

populations=($(cut -f3 --delim="," ${SCRIPTDIR}/3_Results/PredictedPopulations.csv | tail -n +2 | sort | uniq))
for each in ${populations[@]}
do
   grep ${each} ${SCRIPTDIR}/3_Results/PredictedPopulations.csv | cut -f1-2 --delim="," --output-delimiter=" " > ${PROCESSDIR}/${each}Samples.txt
   
   ## 
   if [[ $(wc -l <${PROCESSDIR}/${each}Samples.txt) -ge 1 ]]
   then
      check_relatedness ${PROCESSDIR}/${each}Samples.txt ${each}
   fi
   
done  








