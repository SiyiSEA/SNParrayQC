#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=06aImputationFormatSanger.o
#SBATCH --error=06aImputationFormatSanger.e
#SBATCH --job-name=IFS



## format files for use with Sanger Imputation Server https://imputation.sanger.ac.uk/

## EXECUTION
# sh SNPArray/preprocessing/4_formatForimputation.sh <data path> <population> 
# where 
# <population > is 3 letter code for super population state ALL for no subsetting by population
# <SNP ref file> is an input file of 
# script needs to be executed from <git repo>/array/

## REQUIRES the following variables in config file
# PROCESSDIR, IMPUTEDIR, FILEPREFIX

## REQUIRES the following software
# plink, perl,

## INPUT
#  # binary plink files following prelim QC

## OUTPUT
# vcf files split by chr for upload to michegan imputation server


## IMPORTANT
# before run the imputation, please check the population;
# only for EUR population, the reference panel is HRC;
# anyother population should be 1000G panel.

echo "checking the arguments for config file---------------------------------------------"
datapeth=$1

if [ -z "$1" ]
then
        echo "No argument supplied"
        echo "Please input the paht of the data folder as the first argument"
		    exit 1 # fail
fi

echo "running the PostQCSanger at $datapeth"
source ${datapeth}/config

mv 06aImputationFormatSanger.o ${JOBSDIR}/06aImputationFormatSanger.o
mv 06aImputationFormatSanger.e ${JOBSDIR}/06aImputationFormatSanger.e

touch "$logfile_06a"
exec > >(tee "$logfile_06a") 2>&1
cd ${PROCESSDIR}/FormatImputation || exit 1


echo "checking the arguments--------------------------------------------------"
population=$2

if [ -z "$2" ];
then
        echo "No argument supplied"
        echo "Please input the population argument"
        echo "Population options come from the end of the" $logfile_03b
		exit
else
        if [[ $population == "EUR" ]];
        then
            refFile=${KGG}/../HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab
            mkdir -p InputSangerEUR || exit
		        cd InputSangerEUR || exit

        elif [[ $population == "ALL" ]];
        then
            refFile=${KGG}/1000GP_Phase3_combined.legend
            mkdir -p InputSangerALL || exit
		        cd InputSangerALL || exit
        else
            echo "Cannot process the population" 
            exit
        fi
fi

echo "Start formatting the QC data for Sanger Imputation-----------------------"

module purge
module load VCFtools
module load BCFtools

## use tool to check data prior to upload https://www.well.ox.ac.uk/~wrayner/tools/
# follow the instruction of https://imputation.sanger.ac.uk/?instructions=1#prepareyourdata
## subset samples
if [ $population == "EUR" ]
  then
  ${PLINK}/plink --bfile ${RESULTSDIR}/01/${FILEPREFIX}_QCd_trimmed \
                --keep ${PROCESSDIR}/CheckRelatedness/${population}Samples.txt \
                --maf 0.05 \
                --make-bed \
                --out ${FILEPREFIX}_QCd_${population}
  PCAforPlinkData ${PROCESSDIR}/FormatImputation ${FILEPREFIX}_QCd_${population} 2

elif [ $population == "ALL" ]
  then
    cp ${RESULTSDIR}/02/${FILEPREFIX}_QCd_trimmed.bim ${FILEPREFIX}_QCd_${population}.bim
  	cp ${RESULTSDIR}/02/${FILEPREFIX}_QCd_trimmed.bed ${FILEPREFIX}_QCd_${population}.bed
  	cp ${RESULTSDIR}/02/${FILEPREFIX}_QCd_trimmed.fam ${FILEPREFIX}_QCd_${population}.fam
else
  echo "Please input either EUR or ALL as the first arg!"
  exit 1
fi

echo "Doing liftover -------------------------------------------------------"
# convert the hg17.bim file into hg19.BED file
awk '{print "chr"$1, "\t", $4-1, "\t", $4, "\t", $2}' ${FILEPREFIX}_QCd_${population}.bim > QCd.BED
${LIFTOVER} QCd.BED "${LiftChainHg19}" Mapped.BED unMapped 
mapped_variant=$(wc -l Mapped.BED)
total_variant=$(wc -l QCd.BED)
echo "${mapped_variant}"
echo "${total_variant}"
# echo $(( $mapped_variant*100/$total_variant )) "% of variants have been liftovered successfully!"

# check if you have althernative chr
awk '{print $1}' Mapped.BED | sort -u
awk '{OFS="\t"; print $4, $3}' Mapped.BED > NewPosition.txt

${PLINK}/plink --bfile ${FILEPREFIX}_QCd_${population} \
                --update-map NewPosition.txt \
                --make-bed \
                --out ${FILEPREFIX}_QCd_${population}_hg19


## for HRC check tool need freq file
${PLINK}/plink --bfile ${FILEPREFIX}_QCd_${population}_hg19 \
                --freq \
                --out ${FILEPREFIX}_QCd_${population}_hg19_freq

echo "Running the perl script -------------------------------------------------------"
if [[ $population == "EUR" ]];
then
    refFile=${KGG}/../HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab
    perl ${KINGPATH}/HRC-1000G-check-bim.pl -b ${FILEPREFIX}_QCd_${population}_hg19.bim -f ${FILEPREFIX}_QCd_${population}_hg19_freq.frq -r ${refFile} -h

elif [[ $population == "ALL" ]];
then
    refFile=${KGG}/1000GP_Phase3_combined.legend
    perl ${KINGPATH}/HRC-1000G-check-bim.pl -b ${FILEPREFIX}_QCd_${population}_hg19.bim -f ${FILEPREFIX}_QCd_${population}_hg19_freq.frq -r ${refFile} -g -p ALL
else
    echo "Cannot process the population, Please check the input argument." 
    exit
fi


# sed -i 's=plink=${PLINK}/plink=g' Run-plink.sh
# sed -i '/--real-ref-alleles/d' Run-plink.sh
# sh Run-plink.sh
sed -i 's=plink=${PLINK}/plink=g' Run-plink.sh
head -n 5 Run-plink.sh > Run-plink-Sanger.sh
sh Run-plink-Sanger.sh
${PLINK2} --bfile ${FILEPREFIX}_QCd_${population}_hg19-updated \
          --real-ref-alleles \
          --make-bed \
          --out ${FILEPREFIX}_QCd_${population}_hg19-updated_temp
${PLINK2} --bfile ${FILEPREFIX}_QCd_${population}_hg19-updated_temp \
          --real-ref-alleles \
          --recode vcf \
          --out ${FILEPREFIX}_QCd_${population}_hg19-updated_final
rm TEMP*
rm *temp*

echo "BCFtools -------------------------------------------------------"
bcftools sort ${FILEPREFIX}_QCd_${population}_hg19-updated_final.vcf -Oz -o ${RESULTSDIR}/06a/${FILEPREFIX}_QCd_${population}_hg19_upload.vcf.gz
bcftools index ${RESULTSDIR}/06a/${FILEPREFIX}_QCd_${population}_hg19_upload.vcf.gz


# not works
# sed -i 's=plink=${PLINK}/plink=g' Run-plink.sh
# sh Run-plink.sh
# for file in *.vcf; do vcf-sort ${file} | bgzip -c > ${file}.gz;done
# ls ${FILEPREFIX}_QCd_${population}_hg19-updated-chr*.vcf.gz > list_vsf.txt
# bcftools concat -f list_vsf.txt -Oz -o ${RESULTSDIR}/06a/${FILEPREFIX}_QCd_${population}_hg19.vcf.gz
# bcftools index ${RESULTSDIR}/06a/${FILEPREFIX}_QCd_${population}_hg19.vcf.gz

echo 'done 06a_formatForImputationSanger.sh'