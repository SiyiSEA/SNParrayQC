
## Written by Siyi
## If the --double-id has been applied to convert the vcf file into PLINK format, both family and within-family IDs to be set to the sample ID,
## where the format is FID_IID.
## We cannot seperate the FID_IID simply by '_', because there are lots of '_' in both FID and IID.
## So, this R script only works if FID and IID has too many '_'.

args<-commandArgs(trailingOnly = TRUE)

library(data.table)
filteredfam<-args[1]
rawfam<-args[2]

filter_fam = read.table(filteredfam, header = F)
fam = read.table(rawfam, header = F)

fam$FID_IID = paste0(fam$V1, '_', fam$V2)
#length(setdiff(fam$FID_IID, filter_fam$V1))
m = match(filter_fam[,2], fam$FID_IID)
fam = fam[m,]

write.table(fam[,1:6], filteredfam, col.names = F, row.names = F, quote = F)




