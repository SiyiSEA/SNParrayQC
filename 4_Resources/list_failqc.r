# this script will generate a list of sample fail on the sample-level check

arguments <- commandArgs(T)
result_path <- arguments[1]

setwd(result_path)

if (file.exists(file = "SexCheck.sexcheck") == TRUE){
    sex=read.table("./SexCheck.sexcheck", header = T)
    sex_to_remove = subset(sex, STATUS == "PROBLEM", select = c(FID, IID))
    remove_id = sex_to_remove
}else{
    remove_id = data.frame(FID = NA, IID = NA)
}

if (file.exists(file = "hetGenotypes.het") == TRUE){
    het = read.table("./hetGenotypes.het", header = T)
    het_to_remove = subset(het, F > mean(het$F) + 3*sd(het$F) | F < mean(het$F) - 3*sd(het$F), select = c(FID, IID))
    remove_id = rbind(het_to_remove, remove_id)
}

if (file.exists(file = "missingGenotypes.imiss") == TRUE){
    missing_individual = read.table("./missingGenotypes.imiss", header = T)
    missing_to_remove = subset(missing_individual, F_MISS > 0.01, select = c(FID, IID))
    remove_id = rbind(missing_to_remove, remove_id)
}
remove_id = remove_id[-which(is.na(remove_id)),]
remove_id = unique(remove_id)
message(dim(remove_id)[1], " sample(s) are filed on checking of heterozygosity, missing rate and gender.")

write.table(remove_id, "./fail_mis_het_sex.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)