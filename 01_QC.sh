#!/bin/sh

## DESCRIPTIONS
#This script takes genotype data outputted by genome studio and performs data QC
# including: removal of duplicate samples, update sample ids, perform sex check

## EXECUTION
# sh 01_QC.sh

## SOFTWARES
# king, plink, python, gcta

## REQUIRES the following variables in config file
## INPUT raw_data bfile
## OUTPUT data_QCd PCA 

. ./config

## Sample Duplication QC
cd "${PROCESSDIR}" || exit
mkdir -p QCoutput

# identify duplicates or MZ twins.
"$KINGPATH"/king -b "${RAWDATADIR}"/"${FILEPREFIX}".bed --duplicate --prefix QCoutput/king


if [ -s QCoutput/king.con ]
then 

    cut -f 1,2 QCoutput/king.con > QCoutput/duplicateSamples.tmp
    cut -f 3,4 QCoutput/king.con >> QCoutput/duplicateSamples.tmp
    sort QCoutput/duplicateSamples.tmp | uniq > QCoutput/duplicateSamples.txt
    rm QCoutput/duplicateSamples.tmp

    if [ -s QCoutput/duplicateSamples.txt ]
    then
    
    ## elevated missing data rate
    ${PLINK}/plink --bfile ${RAWDATADIR}/${FILEPREFIX} --missing --out duplicateSamples

	## use python script to identify duplicated with greatest missingness
	python ${SCRIPTDIR}/utilitys/ExcludeDuplicates.py QCoutput/king.con duplicateSamples.imiss QCoutput/dupsToExclude.txt

	## remove duplicates
	${PLINK}/plink --bfile ${RAWDATADIR}/${FILEPREFIX} --remove QCoutput/dupsToExclude.txt --make-bed --out ${FILEPREFIX}_update_1
    fi

else
    cp ${RAWDATADIR}/${FILEPREFIX}.bed ${FILEPREFIX}_update_1.bed
    cp ${RAWDATADIR}/${FILEPREFIX}.bim ${FILEPREFIX}_update_1.bim
    cp ${RAWDATADIR}/${FILEPREFIX}.fam ${FILEPREFIX}_update_1.fam

fi


## Variants Duplication QC - remove variants at the same position (i.e. triallelic)
awk '{if ($1 != 0) print $1":"$4}' ${FILEPREFIX}_update_1.bim > QCoutput/pos.tmp
sort QCoutput/pos.tmp | uniq -d > QCoutput/dupLocs.txt

awk -F ":" '{print $1,$2-1,$2,"set1", "set2"}' QCoutput/dupLocs.txt > QCoutput/positionsExclude.txt
${PLINK}/plink --bfile ${FILEPREFIX}_update_1 --exclude range QCoutput/positionsExclude.txt --make-bed --out ${FILEPREFIX}_update_2

rm QCoutput/pos.tmp
rm QCoutput/dupLocs.txt


## exclude samples which do not have sex predicted
## exclude mismatched samples
## retain samples with missing sex info
${PLINK}/plink --bfile ${FILEPREFIX}_update_2 --mind 0.02 --check-sex --out QCoutput/SexCheck
awk '{if ($4 == 0) print $1,$2 }' QCoutput/SexCheck.sexcheck > QCoutput/sexErrors.txt
awk '{if ($4 != $3 && $3 != 0) print $1,$2 }' QCoutput/SexCheck.sexcheck >> QCoutput/sexErrors.txt
${PLINK}/plink --bfile ${FILEPREFIX}_update_2 --remove QCoutput/sexErrors.txt --make-bed --out ${FILEPREFIX}_update_3



## Sample Heterozygosity and missing genetype call rate QC (only for autosome)
awk '{if ($1 >= 1 && $1 <= 22) print $2}' ${FILEPREFIX}_update_3.bim > QCoutput/autosomalVariants.txt
${PLINK}/plink --bfile ${FILEPREFIX}_update_3 --extract QCoutput/autosomalVariants.txt --missing --out missingSamples

## Smpale related indiciduals (IBS, IBD calculation)
${PLINK}/plink --bfile ${FILEPREFIX}_update_3 --extract QCoutput/autosomalVariants.txt --maf 0.01 --hwe 0.00001 --mind 0.02 --geno 0.05 --indep-pairwise 5000 1000 0.2 --out ld.auto
${PLINK}/plink --bfile ${FILEPREFIX}_update_3 --extract QCoutput/ld.auto.prune.in --het --out QCoutput/roh

## plot missing rate vs F rate, missing rate vs heterozygosity
Rscript /lustre/home/sww208/GoDMC/QC/OrganizedSNParray/4_Resources/Plot_mis_het.R

## exclude anyone with |Fhet| > 0.2, missing rate > 0.02
# based on the mis_het.pdf, the --mind should be 0.02
awk '{if ($6 > 0.2 || $6 < -0.2) print $1,$2}' roh.het > QCoutput/excessHet.txt
${PLINK}/plink --bfile ${FILEPREFIX}_update_3 --remove QCoutput/excessHet.txt --make-bed --out ${FILEPREFIX}_update_4
rm QCoutput/autosomalVariants.txt





## Variants missing call rate QC
# previous QC step for calculating sample missing call rate - .lmiss
# based on the mis_het.pdf, the geno threshold could be 0.05


## Variants Minor allele frequencies
# based on the mis_het.pdf, the maf threshold could be 0.02?
${PLINK}/plink --bfile ${FILEPREFIX}_update_4 --freq --out variant_freq
Rscript /lustre/home/sww208/GoDMC/QC/OrganizedSNParray/4_Resources/Plot_MAF.R

## filter sample and variant missingness, HWE, rare variants and exclude variants with no position
awk '{if ($1 == 0) print $2}' ${FILEPREFIX}_update_4.bim > QCoutput/noLocPos.tmp
${PLINK}/plink --bfile ${FILEPREFIX}_update_4 --exclude QCoutput/noLocPos.tmp --maf 0.001 --hwe 0.00001 --mind 0.02 --geno 0.05 --make-bed --out ${FILEPREFIX}_QCd



## clean up intermediate files but keep log files
rm ${FILEPREFIX}_update_*.b*
rm ${FILEPREFIX}_update_*.fam


## calc PCS within sample only
# LD prune
${PLINK}/plink --bfile ${FILEPREFIX}_QCd --indep 50 5 1.5 --out ${FILEPREFIX}_QCd.ld
${PLINK}/plink --bfile ${FILEPREFIX}_QCd --extract ${FILEPREFIX}_QCd.ld.prune.in --make-bed --out ${FILEPREFIX}_QCd.ld.prune

mkdir -p GCTA

${GCTA}/gcta-1.94.1 --bfile ${FILEPREFIX}_QCd.ld.prune --make-grm-bin --autosome --out GCTA/${FILEPREFIX}_QCd_GCTA
${GCTA}/gcta-1.94.1 --grm GCTA/${FILEPREFIX}_QCd_GCTA --pca --out GCTA/${FILEPREFIX}_QCd.pca

rm ${FILEPREFIX}_QCd.ld.prune*

## extract SNP probes for comparison with DNAm data
${PLINK}/plink --bfile ${FILEPREFIX}_QCd --extract ${EPICREF}/RSprobes.txt --recodeA --out ${FILEPREFIX}_59DNAmSNPs

