---
title: "QCreport"
output: html_notebook
---

# Sample-Level QC vislization
```{r}
# distribution of heterozygosity
setwd("/lustre/home/sww208/QC/SNParrayQC")
het = read.table("2_ProcessedData/QCData/hetGenotypes.het", header = T)
hist(het$F, breaks = 100, col = "lightblue", main = "F = (O-E)/(N-E) of the Heterozygosity", xlab = "F score")
abline(v = mean(het$F), lty = 2)
abline(v = mean(het$F) + 3*sd(het$F), lty = 2, col = "red")
abline(v = mean(het$F) - 3*sd(het$F), lty = 2, col = "red")
```

```{r}
# missing value vs het
missing_individual = read.table("2_ProcessedData/QCData/missingGenotypes.imiss", header = T)
hist(missing_individual$F_MISS, breaks = 100, col = "lightblue", main = "F = (O-E)/(N-E) of the Missingness", xlab = "F score")
abline(v = 0.01, lty = 2, col = "red")
```

```{r}
plot(y=het$F, x=missing_individual$F_MISS, ylab="Heterozygosity rate", xlab="Proportion of missing genotypes")
abline(h = mean(het$F) + 3*sd(het$F), lty = 2, col = "blue")
abline(h = mean(het$F) - 3*sd(het$F), lty = 2, col = "blue")
abline(v = 0.01, lty = 2, col = "red")
```

```{r}
het_to_remove = subset(het, F > mean(het$F) + 3*sd(het$F) | F < mean(het$F) - 3*sd(het$F), select = c(FID, IID))
missing_to_remove = subset(missing_individual, F_MISS > 0.01, select = c(FID, IID))
het_miss_remove = rbind(het_to_remove, missing_to_remove)
het_miss_remove = het_miss_remove[!duplicated(het_miss_remove),]
#write.table(het_miss_remove, "fail_mis_het_qc.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

```{r}
# sex check
sex=read.table("2_ProcessedData/QCData/SexCheck.sexcheck", header = T)
hist(sex$F, breaks = 100, col = "lightblue", main = "F of sex", xlab = "F score")
```
# SNP-Level QC vislization
```{r}
# SNP level
missing_SNP = read.table("2_ProcessedData/QCData/rawMissing.lmiss", header = T)
```

```{r}
hist(log10(missing_SNP$F_MISS), col = "lightblue", main = "Fraction of missing data", xlab = "Fraction of missing genotypes",ylab="Number of SNPs")
abline(v=log10(0.05),lty=2)
```

```{r}
# SNP level
feq_snp = read.table("2_ProcessedData/QCData/rawVariantFreq.frq", header = T)
```

```{r}
hist(feq_snp$MAF, col = "lightblue",breaks=50,ylab="Number of SNPs",xlab="MAF",main="Minor allele frequencies")
abline(v=0.01,lty=2) #1%
abline(v=0.02,lty=2)
```

