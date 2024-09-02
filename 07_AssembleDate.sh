#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=250G
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=07AssembleData.o
#SBATCH --error=07AssembleData.e
#SBATCH --job-name=AssembleData


# this script is to assemble the imputed data
## EXECUTION
# bash ./07_AssembleData.sh

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
touch "$logfile_07"
exec > >(tee "$logfile_07") 2>&1

mv 07AssembleData.o ${JOBSDIR}/07AssembleData.o
mv 07AssembleData.e ${JOBSDIR}/07AssembleData.e

mv 07AssembleData.o ${JOBSDIR}/07AssembleData.o
mv 07AssembleData.e ${JOBSDIR}/07AssembleData.e

echo "Assemble QCd data----------------------------------------------"
cp ${PROCESSDIR}/QCData/${FILEPREFIX}_QCd_trimmed.* ${DATADIR}/3_Results/07/.


echo "Assemble imputed data-------------------------------------------"
assemble_imputed_data () {
	server=$1
	panel=$2

	echo "Assembling imputed data from ${server} for ${panel}"
	cd ${IMPUTEDIR}/ImputationOutput${server}${panel} || exit

	rm -f mergedosefile.txt
	touch mergedosefile.txt

	if [ -s chr_X.zip ]
	then
		chrNum=23
	else
		chrNum=22
	fi

	for i in $(seq 1 $chrNum)
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
	${PLINK2} --bfile data_chr1_dose_temp \
			  --pmerge-list mergedosefile.txt \
			  --make-bed \
			  --out data_dose_${panel}_${server}_temp
	
	# make sure all the chx convert into chr23 across all the files
	awk '{print $2}' data_dose_${panel}_${server}_temp.bim > oldidchrX.txt
	sed 's/X/23/g' oldidchrX.txt > newidchr23.txt
	paste oldidchrX.txt newidchr23.txt > updatechrid.txt
	
	# correct the FID and IID for fam file
	Rscript ${SCRIPTDIR}/4_Resources/correctFIDIID.r data_dose_${panel}_${server}_temp.fam

	${PLINK}/plink --bfile data_dose_${panel}_${server}_temp \
				--make-bed \
				--update-ids updateFIDIID.txt \
				--update-name updatechrid.txt 2 1 \
				--out data_dose_${panel}_${server}

	# remove reducdecy data
	rm data_dose_*temp* mergedosefile.txt
	rm oldidchrX.txt newidchr23.txt updatechrid.txt updateFIDIID.txt
	cp data_dose_${panel}_${server}.*  ${DATADIR}/3_Results/07/.
}

# pass
# assemble_imputed_data Michigan HRC
assemble_imputed_data Sanger HRC
# assemble_imputed_data Michigan 1000G
# assemble_imputed_data Sanger 1000G

echo "Assemble filtered data------------------------------------------"
assemble_imputed_postQC_data () {
	server=$1
	panel=$2
	echo "Assembling Filtereed Data from ${server} for ${panel}"
	cd ${IMPUTEDIR}/ImputationOutput${server}${panel} || exit
	pwd
	cp data_filtered_${server}.bim  data_filtered_${panel}_${server}.bim 
	cp data_filtered_${server}.bed  data_filtered_${panel}_${server}.bed
	cp data_filtered_${server}.fam  data_filtered_${panel}_${server}.fam
	cp data_filtered_${server}.info  data_filtered_${panel}_${server}.info
	mv data_filtered_${panel}_${server}* ${DATADIR}/3_Results/07/.
}

# pass
# assemble_imputed_postQC_data Michigan HRC
assemble_imputed_postQC_data Sanger HRC
# assemble_imputed_postQC_data Michigan 1000G
# assemble_imputed_postQC_data Sanger 1000G

echo "Plotting filtered data------------------------------------------"
plot_imputed_postQC_data () {
	server_1=Michigan
	server_2=Sanger
	panel=$1

	if [ -s ${DATADIR}/3_Results/07/data_filtered_${panel}_${server_1}.info ]
	then
		if [ -s ${DATADIR}/3_Results/07/data_filtered_${panel}_${server_2}.info ]
		then
			echo "Plotting the density plots for ${server_1}${panel} and ${server_2}${panel}"
			Rscript ${SCRIPTDIR}/4_Resources/plot_imputation_quality.r \
				"${DATADIR}/3_Results/07/" \
				"./data_filtered_${panel}_${server_1}.info" \
				"./data_filtered_${panel}_${server_2}.info" \
				${panel}
		else
			echo "Plotting the density plots for ${server_1}${panel}"
			Rscript ${SCRIPTDIR}/4_Resources/plot_imputation_quality.r \
				"${DATADIR}/3_Results/07/" \
				"./data_filtered_${panel}_${server_1}.info" \
				"NULL" \
				${panel}
		fi
	else
		if [ -s ${DATADIR}/3_Results/07/data_filtered_${panel}_${server_2}.info ]
		then
			echo "Plotting the density plots for ${server_1}${panel} and ${server_2}${panel}"
			Rscript ${SCRIPTDIR}/4_Resources/plot_imputation_quality.r \
				"${DATADIR}/3_Results/07/" \
				"NULL" \
				"./data_filtered_${panel}_${server_2}.info" \
				${panel}
		else
			echo "Error: There is no valid file of info for plotting."
		fi
	fi

}

plot_imputed_postQC_data HRC 
# plot_imputed_postQC_data 1000G


echo "PCA and count variants------------------------------------------"
# PCA and count the variants
PCA_count_variant () {
	data_path=$1
	data_prefix=$2
	PCApath=${DATADIR}/3_Results/PCAVariants
	
	echo "Plot PCA and count the outliers for ${data_prefix}"
	# PCA
	${PLINK}/plink --bfile ${data_path}/${data_prefix} --indep 50 5 1.5 --out ${PCApath}/${data_prefix}.ld
	${PLINK}/plink --bfile ${data_path}/${data_prefix} --extract ${PCApath}/${data_prefix}.ld.prune.in --make-bed --out ${PCApath}/${data_prefix}.ld.prune

	${GCTA}/gcta-1.94.1 --bfile ${PCApath}/${data_prefix}.ld.prune --make-grm-bin --autosome --out ${PCApath}/${data_prefix}.imqc
	${GCTA}/gcta-1.94.1 --grm ${PCApath}/${data_prefix}.imqc --pca --thread-num 4 --out ${PCApath}/${data_prefix}.imqc.pca

	# plot PCs to identify outliers
	Rscript ${SCRIPTDIR}/4_Resources/plotPCs.r ${PCApath}/${data_prefix}.imqc.pca.eigenvec 3

	wc -l ${data_path}/${data_prefix}.bim >> ${PCApath}/VariantsCount.txt

	for nPC in PC1 PC2 PC3 PC4 PC5 PC6
	do
		num_outlier=$(awk -v npc=${nPC} '$4==npc' ${PCApath}/${data_prefix}_OutliersFromPC_3SDfromMean.txt | wc -l)
		echo  ${nPC} ${num_outlier} ${PCApath}/${data_prefix}_OutliersFromPC_3SDfromMean.txt >> ${PCApath}/VariantsCount.txt
	done

	rm ${PCApath}/${data_prefix}*.ld* 
	rm ${PCApath}/*.imqc* 

	awk 'NR>1' ${PCApath}/${data_prefix}_OutliersFromPC_3SDfromMean.txt > ${PCApath}/${data_prefix}_OutliersFromPC_3SDfromMean_temp.txt
	sed -e 'ROW FID IID TYPE nPC' ${PCApath}/${data_prefix}_OutliersFromPC_3SDfromMean_temp.txt
}


cd ${DATADIR}/3_Results/PCAVariants || exit
rm -f VariantsCount.txt
touch VariantsCount.txt

# takes a while
PCA_count_variant ${DATADIR}/3_Results/07 data_dose_HRC_Sanger
PCA_count_variant ${DATADIR}/3_Results/07 data_filtered_HRC_Sanger


rm ${RESULTSDIR}/PCAVariants/${FILEPREFIX}*.eigenv*
rm ${RESULTSDIR}/PCAVariants/${FILEPREFIX}*.impc.grm*