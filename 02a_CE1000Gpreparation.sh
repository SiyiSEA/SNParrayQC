#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=/lustre/home/sww208/QC/SNParrayQC/5_JobReports/02a1000Gpreparation.o
#SBATCH --error=/lustre/home/sww208/QC/SNParrayQC/5_JobReports/02a1000Gpreparation.e
#SBATCH --job-name=02a

# This scritp is to generating 1000G based on hg19;
# if you do not have the 1000G in hg38, you can run this script;
source ./config
touch "$logfile_02"
source ${RESOURCEDIR}/PCAforPlinkData.sh
exec > >(tee "$logfile_02") 2>&1
cd ${RESOURCEDIR}/1000G || exit

#================================================================================================
# The first 1000G comes from https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3

pgen=https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1
pvar=https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst?dl=1
sample=https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?rlkey=0yyifzj2fb863ddbmsv4jkeq6&dl=1

echo "download the pgen file-----------------------------------------------------------------"
wget $pgen
mv 'all_phase3.pgen.zst?dl=1' all_phase3.pgen.zst
${PLINK2} --zst-decompress all_phase3.pgen.zst > all_phase3.pgen

echo "download the pvar file-----------------------------------------------------------------"
wget $pvar
mv 'all_phase3.pvar.zst?dl=1' all_phase3.pvar.zst
${PLINK2} --zst-decompress all_phase3.pvar.zst > all_phase3.pvar

echo "download the psam file-----------------------------------------------------------------"
wget $sample
mv 'phase3_corrected.psam?rlkey=0yyifzj2fb863ddbmsv4jkeq6' all_phase3.psam

echo "convert into PLINK binary files--------------------------------------------------------"
${PLINK2} --pfile all_phase3 vzs \
        --max-alleles 2 \
        --make-bed \
        --autosome \
        --out all_phase3_1000G

# check the chrom
awk '{print $1}' all_phase3_1000G.bim | sort | uniq

awk '{print $1":"$4_$5_$6, $2}' all_phase3_1000G.bim > 1000GVariantIDConvert.txt

#  Prune variants
${PLINK}/plink --bfile all_phase3_1000G \
                --maf 0.01 \
                --geno 0.01 \
                --make-bed \
                --out all_phase3_1000G_update_2

# rename the vairants id
${PLINK2} --bfile all_phase3_1000G_update_2 \
                --set-all-var-ids @:#_\$1_\$2 \
                --make-bed \
                --new-id-max-allele-len 662 \
                --out all_phase3_1000G_update_3

echo "Please set the Ref1000G as ${RESOURCEDIR}/1000G/all_phase3_1000G_update_3 in the config file"
echo "Do the liftover for your data to hg19 --> set your liftoverChain to the matched version."