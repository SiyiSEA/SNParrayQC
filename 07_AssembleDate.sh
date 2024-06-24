#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=2 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=250G
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=/lustre/home/sww208/QC/SNParrayQC/5_JobReports/AssembleData.o
#SBATCH --error=/lustre/home/sww208/QC/SNParrayQC/5_JobReports/AssembleData.e
#SBATCH --job-name=AssembleData


# this script is to assemble the imputed data
## EXECUTION
# sh ./06_AssembleData.sh

## REQUIRES the following variables in config file
# ${IMPUTEDIR}, ${RAWDATADIR}/${FILEPREFIX}, ${SCRIPTDIR}

## REQUIRES the following software
# plink1.9, plink2, GCTA, R

## INPUT
# raw data
# QCd data
# Imputed data on Michigan and Sanger
# ImputedPost data on Michigan and Sanger

## OUTPUT
# four folders in the 3_Results
# three of them are for QCdata, ImputedData and ImputedPostData
# one is PCAVariants for sorting PCA plots, outliers variants and number of variants



source ./config
echo 'runing 07_AssembleData.sh'
module load R/4.2.1-foss-2022a

echo "Assemble QCd data----------------------------------------------"
# assemble_QC_data 
# mkdir -p ${SCRIPTDIR}/3_Results/QCdData
# cd ${PROCESSDIR} || exit
# cp ${FILEPREFIX}_QCd.b* ${FILEPREFIX}_QCd.fam ${SCRIPTDIR}/3_Results/QCdData/.


echo "Assemble imputed data-------------------------------------------"
assemble_imputed_data () {
	server=$1
	panel=$2

	echo "Assembling imputed data from ${server} for ${panel}"
	cd ${SCRIPTDIR}/3_Results/ || exit
	mkdir -p ${SCRIPTDIR}/3_Results/ImputedData${server}${panel}
	cd ${IMPUTEDIR}/ImputationOutput${server}${panel} || exit

	rm -f mergedosefile.txt
	touch mergedosefile.txt

	for i in {1..23}
	do 
	
		if [[ $panel = "1000G" && $i -eq 23 ]]
		then
			${PLINK2} --pfile data_chr${i}_dose \
				  --rm-dup exclude-mismatch \
				  --make-pgen \
				  --merge-par \
				  --out data_chr${i}_dose_temp
		else
			${PLINK2} --pfile data_chr${i}_dose \
				  --rm-dup exclude-mismatch \
				  --make-pgen \
				  --out data_chr${i}_dose_temp
		fi
		${PLINK2} --pfile data_chr${i}_dose_temp \
			  	  --make-bed \
			  	  --out data_chr${i}_dose_temp

		echo "data_chr${i}_dose_temp" >>"mergedosefile.txt"
	done

	sed -i '1d' mergedosefile.txt
	${PLINK2} --bfile data_chr1_dose_temp --pmerge-list mergedosefile.txt --make-bed --out data_dose_${panel}_${server}_temp
	
	# make sure all the chx convert into chr23 across all the files
	awk '{print $2}' data_dose_${panel}_${server}_temp.bim > oldidchrX.txt
	sed 's/X/23/g' oldidchrX.txt > newidchr23.txt
	paste oldidchrX.txt newidchr23.txt > updatechrid.txt
	
	${PLINK}/plink --bfile data_dose_${panel}_${server}_temp --make-bed --update-name updatechrid.txt 2 1 --out data_dose_${panel}_${server}

	# correct the FID and IID for fam file
	cp data_dose_${panel}_${server}.fam data_dose_${panel}_${server}.fam.orig
	awk '{print $2,$2,$3,$4,$5,$6}' < data_dose_${panel}_${server}.fam.orig > data_dose_${panel}_${server}.fam

	Rscript ${DATADIR}/4_Resources/correctFIDIID.R \
		data_dose_${panel}_${server}.fam \
		${RAWDATADIR}/${FILEPREFIX}.fam

	# remove reducdecy data
	rm data_dose_*temp* data_dose_${panel}_${server}.fam.orig mergedosefile.txt
	rm oldidchrX.txt newidchr23.txt updatechrid.txt
	cp data_dose_${panel}_${server}.b* data_dose_${panel}_${server}.fam  ${SCRIPTDIR}/3_Results/ImputedData${server}${panel}/.
}

#pass
# assemble_imputed_data Michigan HRC
# assemble_imputed_data Sanger HRC
# assemble_imputed_data Michigan 1000G
# assemble_imputed_data Sanger 1000G

echo "Assemble filtered data------------------------------------------"
assemble_imputed_postQC_data () {
	server=$1
	panel=$2
	echo "Assembling Filtereed Data from ${server} for ${panel}"
	mkdir -p ${SCRIPTDIR}/3_Results/FilterData${server}${panel}
	cd ${IMPUTEDIR}/ImputationOutput${server}${panel} || exit
	pwd
	cp data_filtered_${server}.bim  data_filtered_${panel}_${server}.bim 
	cp data_filtered_${server}.bed  data_filtered_${panel}_${server}.bed
	cp data_filtered_${server}.fam  data_filtered_${panel}_${server}.fam
	cp data_filtered_${server}.info  data_filtered_${panel}_${server}.info
	mv data_filtered_${panel}_${server}* ${SCRIPTDIR}/3_Results/FilterData${server}${panel}/.
}

# pass
# assemble_imputed_postQC_data Michigan HRC
# assemble_imputed_postQC_data Sanger HRC
# assemble_imputed_postQC_data Michigan 1000G
# assemble_imputed_postQC_data Sanger 1000G

echo "Plotting filtered data------------------------------------------"

plot_imputed_postQC_data () {
	server_1=Michigan
	server_2=Sanger
	panel=$1

	echo "Plotting the density plots for ${server_1}${panel} and ${server_2}${panel}"
	Rscript ${DATADIR}/4_Resources/plot_imputation_quality.R \
		"${DATADIR}/3_Results/" \
		"./FilterData${server_1}${panel}/data_filtered_${panel}_${server_1}.info" \
		"./FilterData${server_2}${panel}/data_filtered_${panel}_${server_2}.info" \
		${panel}

}

# plot_imputed_postQC_data HRC 
# plot_imputed_postQC_data 1000G


echo "PCA and count variants------------------------------------------"
# PCA and count the variants
PCA_count_variant () {
	data_path=$1
	data_prefix=$2
	PCApath=${SCRIPTDIR}/3_Results/PCAVariants
	
	echo "Plot PCA and count the outliers for ${data_prefix}"
	# PCA
	${PLINK}/plink --bfile ${data_path}/${data_prefix} --indep 50 5 1.5 --out ${PCApath}/${data_prefix}.ld
	${PLINK}/plink --bfile ${data_path}/${data_prefix} --extract ${PCApath}/${data_prefix}.ld.prune.in --make-bed --out ${PCApath}/${data_prefix}.ld.prune

	${GCTA}/gcta-1.94.1 --bfile ${PCApath}/${data_prefix}.ld.prune --make-grm-bin --autosome --out ${PCApath}/${data_prefix}.imqc
	${GCTA}/gcta-1.94.1 --grm ${PCApath}/${data_prefix}.imqc --pca --out ${PCApath}/${data_prefix}.imqc.pca

	# plot PCs to identify outliers
	Rscript ${SCRIPTDIR}/4_Resources/plotPCs.r ${PCApath}/${data_prefix}.imqc.pca.eigenvec 3

	wc -l ${data_path}/${data_prefix}.bim >> ${PCApath}/VariantsCount.txt

	for nPC in PC1 PC2 PC3 PC4 PC5 PC6
	do
		num_outlier=$(awk -v npc=${nPC} '$5==npc' ${PCApath}/${data_prefix}_OutliersFromPC_3SDfromMean.txt | wc -l)
		echo  ${nPC} ${num_outlier} ${PCApath}/${data_prefix}_OutliersFromPC_3SDfromMean.txt >> ${PCApath}/VariantsCount.txt
	done

	rm ${PCApath}/${data_prefix}*.ld* 
	rm ${PCApath}/*.imqc* 

	awk 'NR>1' ${PCApath}/${data_prefix}_OutliersFromPC_3SDfromMean.txt > ${PCApath}/${data_prefix}_OutliersFromPC_3SDfromMean_temp.txt
	sed -e 'ROW FID IID TYPE nPC' ${PCApath}/${data_prefix}_OutliersFromPC_3SDfromMean_temp.txt
}


# cd ${SCRIPTDIR}/3_Results || exit
# mkdir -p PCAVariants
# cd PCAVariants || exit
# rm -f VariantsCount.txt
# touch VariantsCount.txt

# # takes an hour
# PCA_count_variant ${RAWDATADIR} ${FILEPREFIX}
# PCA_count_variant ${SCRIPTDIR}/3_Results/QCdData ${FILEPREFIX}_QCd
# PCA_count_variant ${SCRIPTDIR}/3_Results/ImputedDataSangerHRC data_dose_HRC_Sanger
# PCA_count_variant ${SCRIPTDIR}/3_Results/ImputedDataMichiganHRC data_dose_HRC_Michigan
# PCA_count_variant ${SCRIPTDIR}/3_Results/ImputedDataSanger1000G data_dose_1000G_Sanger
# PCA_count_variant ${SCRIPTDIR}/3_Results/ImputedDataMichigan1000G data_dose_1000G_Michigan
# PCA_count_variant ${SCRIPTDIR}/3_Results/FilterDataMichigan1000G data_filtered_1000G_Michigan
# PCA_count_variant ${SCRIPTDIR}/3_Results/FilterDataSanger1000G data_filtered_1000G_Sanger
# PCA_count_variant ${SCRIPTDIR}/3_Results/FilterDataMichiganHRC data_filtered_HRC_Michigan
# PCA_count_variant ${SCRIPTDIR}/3_Results/FilterDataSangerHRC data_filtered_HRC_Sanger

echo "Relationship between the PCs between QCd data and PostQC data-----------------------------"
