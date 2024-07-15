#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=/lustre/home/sww208/QC/SNParrayQC/ImputationFormatMichigen.o
#SBATCH --error=/lustre/home/sww208/QC/SNParrayQC/ImputationFormatMichigen.e
#SBATCH --job-name=IFM


## format files for use with Michegan Imputation Server

## EXECUTION
# sbatch 05a_formatForimputation.sh <population> <SNP ref file>
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


source ./config
touch "$logfile_05a"
exec > >(tee "$logfile_05a") 2>&1
cd ${PROCESSDIR}/FormatImputation || exit

echo "checking the arguments--------------------------------------------------"
population=$1

if [ -z "$1" ];
then
        echo "No argument supplied"
        echo "Please input the population argument"
        echo "Population options come from the end of the" $logfile_03b
		exit
else
        if [[ $population == "EUR" ]];
        then
            refFile=${KGG}/../HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab
            mkdir -p InputMichiganEUR || exit
		    cd InputMichiganEUR || exit

        elif [[ $population == "ALL" ]];
        then
            refFile=${KGG}/1000GP_Phase3_combined.legend
            mkdir -p InputMichiganALL || exit
		    cd InputMichiganALL || exit
        else
            echo "Cannot process the population" 
            exit
        fi
fi

echo "Start formatting the QC data for Michigan Imputation-----------------------"

module purge
module load VCFtools

## use tool to check data prior to upload https://www.well.ox.ac.uk/~wrayner/tools/
# folow the instruction from https://imputationserver.readthedocs.io/en/latest/prepare-your-data/

## subset samples
if [ $population == "EUR" ]
  then
  ${PLINK}/plink --bfile${RESULTSDIR}/02/${FILEPREFIX}_QCd_trimmed \
                --keep ${PROCESSDIR}/${population}Samples.txt \
                --maf 0.05 \
                --make-bed \
                --out ${FILEPREFIX}_QCd_${population}
elif [ $population == "ALL" ]
  then
    cp ${RESULTSDIR}/02/${FILEPREFIX}_QCd_trimmed.bim ${FILEPREFIX}_QCd_${population}.bim
  	cp ${RESULTSDIR}/02/${FILEPREFIX}_QCd_trimmed.bed ${FILEPREFIX}_QCd_${population}.bed
  	cp ${RESULTSDIR}/02/${FILEPREFIX}_QCd_trimmed.fam ${FILEPREFIX}_QCd_${population}.fam
else
  echo "Please input either EUR or ALL as the first arg!"
  exit 1
fi

# convert the hg17.bim file into hg19.BED file
awk '{print "chr"$1, "\t", $4-1, "\t", $4, "\t", $2}' ${FILEPREFIX}_QCd_${population} > QCd.BED
${LIFTOVER} QCd.BED "${LiftChainHg19}" Mapped.BED unMapped 
mapped_variant=$(wc -l Mapped.BED)
total_variant=$(wc -l QCd.BED)
echo $(( $mapped_variant*100/$total_variant )) "% of variants have been liftovered successfully!"

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

sed -i 's=plink=${PLINK}/plink=g' Run-plink.sh
sh Run-plink.sh

for file in *.vcf; do vcf-sort ${file} | bgzip -c > ${RESULTSDIR}/05a/${file}.gz;done
rm *.vcf
rm ${FILEPREFIX}_QCd*.*[^gz]
echo 'done 05a_formatForImputation.sh'
