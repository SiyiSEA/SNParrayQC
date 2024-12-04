library(data.table)

print("Reading data")

info = fread("EXTEND.all.info", header = F, data.table=F)

head(info)

fre = fread("EXTEND_imputed_928_sex_updateID_freq.frq", header = T, data.table = F)

head(fre)

print("Matching info and frq")

info_filter = subset(info, info$V1 %in% fre$SNP)

#info_filter_no_dup = info_filter[!duplicated(info_filter$V1), ]

head(info_filter)

fre_filter = subset(fre, fre$SNP %in% info_filter$V1, select =c("SNP", "MAF") )

info_filter2 = merge(fre_filter, info_filter, by.x = "SNP", by.y = "V1")

head(info_filter2)

colnames(info_filter2) = c("SNP", "MAF", "R2")


print("Saving info file as EXTEND.info and saving filtered variants IDs.")

write.table(info_filter2, "EXTEND.info", quote = F, sep = " ", row.names = F)

write.table(fre_filter$SNP, file = "EXTEND_variants_keep.txt", quote = F, sep = " ", row.names = F, col.names = F)
