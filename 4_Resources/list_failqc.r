# this script will generate a list of sample fail on the sample-level check

arguments <- commandArgs(T)
result_path <- arguments[1]

setwd(result_path)

het = read.table("./hetGenotypes.het", header = T)
missing_individual = read.table("./missingGenotypes.imiss", header = T)
sex=read.table("./SexCheck.sexcheck", header = T)

het_to_remove = subset(het, F > mean(het$F) + 3*sd(het$F) | F < mean(het$F) - 3*sd(het$F), select = c(FID, IID))
missing_to_remove = subset(missing_individual, F_MISS > 0.01, select = c(FID, IID))
sex_to_remove = subset(sex, STATUS == "PROBLEM", select = c(FID, IID))

remove_id = rbind(het_to_remove, missing_to_remove, sex_to_remove)
remove_id = unique(remove_id)
message(dim(remove_id)[1], " sample(s) are filed on checking of heterozygosity, missing rate and gender.")

write.table(remove_id, "./fail_mis_het_sex.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)