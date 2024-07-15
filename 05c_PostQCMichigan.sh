#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=/lustre/home/sww208/QC/SNParrayQC/5_JobReports/PostQCMichigan.o
#SBATCH --error=/lustre/home/sww208/QC/SNParrayQC/5_JobReports/PostQCMichigan.e
#SBATCH --job-name=PostQCMichigan


## output files for use with Michegan Imputation Server
# this script is for deal with the output file from Michegan Imputation Server
# So, please run the imputation on Michegan first and download the zip file for each chr into the ${IMPUTEDIR}/chrzip folder
# keep the passcode that Michegan server send to you!

## EXECUTION
# sbatch ./05b_PostQCMichigan.sh <reference panel> <passcode>

## REQUIRES the following variables in config file
# ${IMPUTEDIR}, ${RAWDATADIR}/${FILEPREFIX}, ${SCRIPTDIR}

## REQUIRES the following software
# plink1.9, plink2, GCTA, R

## INPUT
# chr.zip under the  ${IMPUTEDIR}/chrzip folder

## OUTPUT
# data_filtered_Michigan.bim, data_filtered_Michigan.bed, data_filtered_Michigan.fam data_filtered_Michigan.info

module purge
module load R/4.2.1-foss-2022a
# module load p7zip
source ./config

# function for unzip all the zip file
un7zip () {
	passcode=$1
	for chr in $(seq 1 22); do
		7z x -p"${passcode}" ./chr_$chr.zip;
		gzip -d chr${chr}.info.gz
	done
	7z x -p"${passcode}" ./chr_X.zip
	gzip -d chrX.info.gz
	mv chrX.dose.vcf.gz chr23.dose.vcf.gz
	mv chrX.info chr23.info
}

echo "checking the arguments--------------------------------------------------"
cd ${IMPUTEDIR} || exit

panel=$1

if [ -z "$1" ];
then
        echo "No argument supplied"
        echo "Please input the argument either 1000G nor HRC"
		exit
else
		mkdir -p ImputationOutputMichigan${panel} || exit
		cd ImputationOutputMichigan${panel} || exit
fi

# unzip the result
if [ -z "$2" ];
then
        echo "No passcode supplied"
        echo "Skip the setp for unziping the chr.zip"
else
		passcode=$2
		echo "The passcode is" ${passcode}
        echo "chr.zip files are being unzipped"
		un7zip ${passcode}
fi

echo 'runing PostQC for for imputed data from Michigan------------------------'

# chr23 will fail on --make-pgen if there is no sex info
awk '{print $1,$2,$5}' ${RAWDATADIR}/${FILEPREFIX}.fam > sex.info

for i in {1..23}
do
	# convert the vcf.gz into pgen formats
	if [[ $i -eq 23  && $panel == "1000G" ]]; 
	then
	# chrX in present in the input files, needs sex.info
	# Human chrX pseudoautosomal variant(s) appear in chrX on 1000G, needs split-par or merge-par
		${PLINK2} --vcf chr${i}.dose.vcf.gz \
				--make-pgen \
				--update-sex sex.info \
				--out data_chr${i}_dose \
				--split-par 'hg19'

	elif [[ $i -eq 23 && $panel == "HRC" ]];
	then
		${PLINK2} --vcf chr${i}.dose.vcf.gz \
				--make-pgen \
				--update-sex sex.info \
				--out data_chr${i}_dose
	else
		${PLINK2} --vcf chr${i}.dose.vcf.gz \
				  --make-pgen \
				  --out data_chr${i}_dose
	fi

	# Keep SNPs with MAF>0.01, R2>0.8. hwe 1e-5
	${PLINK2} --pfile data_chr${i}_dose \
			  --extract-if-info "R2>0.8" \
			  --maf 0.01 \
			  --hwe 1e-5 \
			  --make-pgen \
			  --out data_chr${i}_filtered_temp2\
			  --keep-allele-order

	if [[ $i -eq 23 && $panel == "1000G" ]]; 
	then
		# Remove duplicate variants
		${PLINK2} --pfile data_chr${i}_filtered_temp2 \
				--rm-dup exclude-mismatch \
				--make-pgen \
				--out data_chr${i}_filtered \
				--keep-allele-order
	else 
		# Remove duplicate variants
		${PLINK2} --pfile data_chr${i}_filtered_temp2 \
				--rm-dup exclude-mismatch \
				--make-pgen \
				--out data_chr${i}_filtered \
				--keep-allele-order
	fi

	# Convert `pvar/pgen/psam` to `bim/bed/fam` format
	${PLINK2} --pfile data_chr${i}_filtered \
			  --make-bed \
			  --out data_chr${i}_filtered

	#Reformat pvar file to id maf r2 
	awk -v x='#CHROM' '$0~x,EOF''{print $3,$7}' < data_chr${i}_filtered.pvar | perl -pe 's/;/ /g'|awk '{print $1,$3,$4}' |perl -pe 's/ID/ID MAF R2/g' |perl -pe 's/MAF=//g'| perl -pe 's/R2=//g' > data_chr${i}_filtered_temp4.info
	sed 's/X/23/g' data_chr${i}_filtered_temp4.info > data_chr${i}_filtered.info
	
done

# Merge them into one dataset
echo "Merge--------------------------------------------------"
rm -f mergefile.txt
touch mergefile.txt
for i in {2..23}
do 
	if [ -f "data_chr${i}_filtered.bed" ]; then
		echo "data_chr${i}_filtered" >> mergefile.txt
	fi
done

${PLINK2} --bfile data_chr1_filtered \
		  --pmerge-list mergefile.txt \
		  --make-bed \
		  --out data_filtered_Michigan_temp

awk '{print $2}' data_filtered_Michigan_temp.bim > oldidchrX.txt
sed 's/X/23/g' oldidchrX.txt > newidchr23.txt
paste oldidchrX.txt newidchr23.txt > updatechrid.txt

# make sure the all the chrX convert to chr23
${PLINK}/plink --bfile data_filtered_Michigan_temp \
			  --make-bed \
			  --update-name updatechrid.txt 2 1 \
			  --out data_filtered_Michigan

# The FIDs in the fam file are set to 0.
cp data_filtered_Michigan.fam data_filtered_Michigan.fam.orig
awk '{print $2,$2,$3,$4,$5,$6}' < data_filtered_Michigan.fam.orig > data_filtered_Michigan.fam

# correct the FID and IID for fam file
Rscript ${DATADIR}/4_Resources/correctFIDIID.R \
	data_filtered_Michigan.fam \
	${RAWDATADIR}/${FILEPREFIX}.fam

# Combine info files into a single file
cp data_chr1_filtered.info data_filtered_Michigan.info

for i in {2..23}
do
	if [ -f "data_chr${i}_filtered.bed" ]; then
		awk ' NR>1 {print $0}' < data_chr${i}_filtered.info | cat >> data_filtered_Michigan.info
	fi
done

# clean up redundant files
rm data_chr*_filtered_temp* data_chr*_filtered.p*
rm data_filtered_Michigan.fam.orig
rm data_filtered_Michigan_temp* oldidchrX.txt newidchr23.txt updatechrid.txt
rm data_chr*_filtered.info

# check whether the number of variants is the same
echo "The number of variants in bim file is" $(wc -l data_filtered_Michigan.bim)
echo "The number of variants in info file is" $(wc -l data_filtered_Michigan.info)
