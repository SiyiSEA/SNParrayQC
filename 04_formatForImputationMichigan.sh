#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=/lustre/home/sww208/QC/OrganizedSNParray/ImputationFormatMichigen.o
#SBATCH --error=/lustre/home/sww208/QC/OrganizedSNParray/ImputationFormatMichigen.e
#SBATCH --job-name=IFM



## format files for use with Michegan Imputation Server

## EXECUTION
# sh SNPArray/preprocessing/4_formatForimputation.sh <population> <SNP ref file>
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

echo 'runing 4_formatForImputation.sh'

module purge
module load VCFtools
# sh $homedir/SNPArray/preprocessing/4_formatForImputation.sh ALL ${KGG}/1000GP_Phase3_combined.legend
# sh $homedir/SNPArray/preprocessing/4_formatForImputation.sh EUR ${KGG}/../HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab



population=EUR
refFile=${KGG}/../HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab

cd ${IMPUTEDIR}/

mkdir -p ImputationInput

cd ImputationInput

mkdir -p ${population}



## use tool to check data prior to upload https://www.well.ox.ac.uk/~wrayner/tools/
## for All use 1000G
cd ${population}/

## subset samples
${PLINK}/plink --bfile ${PROCESSDIR}/${FILEPREFIX}_QCd --keep ${PROCESSDIR}/${population}Samples.txt --maf 0.05 --out ${FILEPREFIX}_QCd_${population} --make-bed


## liftover to hg19 for imputation
# convert the hg18.bim file into hg18.BED file
cd ${PROCESSDIR} || exit
awk '{print "chr"$1, "\t", $4-1, "\t", $4, "\t", $2}' ${FILEPREFIX}_QCd.bim > hg17.BED
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg17/liftOver/hg17ToHg19.over.chain.gz
${LIFTOVER} hg17.BED hg17ToHg19.over.chain.gz Hg19.BED unMappedHg19

# check if you have althernative chr
awk '{print $1}' Hg19.BED | sort -u
awk '{OFS="\t"; print substr($1,4), $4, 0, $3}' Hg19.BED > Hg19.bim
awk '{OFS="\t"; print $4, $3}' Hg19.BED > Hg19build.txt


cd ${IMPUTEDIR}/ImputationInput/${population}/
${PLINK}/plink --bfile ${FILEPREFIX}_QCd_${population} --update-map ${PROCESSDIR}/Hg19build.txt --make-bed --out ${FILEPREFIX}_QCd_${population}_hg19

## for HRC check tool need freq file
${PLINK}/plink --bfile ${FILEPREFIX}_QCd_${population}_hg19 --freq --out ${FILEPREFIX}_QCd_hg19_freq
perl ${KINGPATH}/HRC-1000G-check-bim.pl -b ${FILEPREFIX}_QCd_${population}_hg19.bim -f ${FILEPREFIX}_QCd_hg19_freq.frq -r ${refFile} -g --hrc


#sed -i 's=plink=${PLINK}/plink=g' Run-plink.sh
#sh Run-plink.sh

#for file in *.vcf; do vcf-sort ${file} | bgzip -c > ${file}.gz;done
#rm *.vcf
#rm ${FILEPREFIX}_QCd*.*[^gz]
echo 'done 4_formatForImputation.sh'
