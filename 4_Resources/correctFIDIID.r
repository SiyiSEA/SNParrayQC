library(stringr)
## Written by Siyi
## Sanger combined the FID and IID with '_'
## Treats the FID_IID both as FID and IID, which is not match to the raw fam file
## this script will generate a txt file for update the FID and IID

args<-commandArgs(trailingOnly = TRUE)

library(data.table)
filteredfam<-args[1]

filter_fam = read.table(filteredfam, header = F)
filter_fam[c("newFID", "newIID")] = str_split_fixed(filter_fam$V2, '_', 2)

write.table(filter_fam[,c("V1", "V2", "newFID", "newIID")], file = "updateFIDIID.txt", col.names = F, row.names = F, quote = F)


