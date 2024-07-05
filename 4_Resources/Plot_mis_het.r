setwd('//lustre/home/sww208/GoDMC/QC/OrganizedSNParray/2_ProcessedData')


# arguments <- commandArgs(T)
# formattedData <- arguments[1]
# clocksDir <- arguments[2]
# agePred <- arguments[3]

het = read.table("./QCoutput/roh.het", header = T)
mis = read.table("./QCoutput/missingSamples.imiss", header = T)
lmis = read.table("./QCoutput/missingSamples.lmiss", header = T)
mishet=data.frame(FID=het$FID, IID=het$IID,
                  het.rate=(het$N.NM - het$O.HOM)/het$N.NM, mis.rate=mis$F_MISS)
Fshet = data.frame(FID=het$FID, IID=het$IID, F=het$F, mis.rate=mis$F_MISS)

#open .pdf file
pdf("QCoutput/mis_het.pdf")
#plot het.rate vs. mis.rate
plot(y=mishet$het.rate , x=mishet$mis.rate, ylab="Individual Heterozygosity rate", xlab="Proportion of missing genotypes")
abline(v=0.01, lty=2)
abline(v=0.02, lty=2)
abline(h=0.305, lty=2)


plot(y=Fshet$F, x=Fshet$mis.rate, ylab="Individual F rate", xlab="Proportion of missing genotype")
abline(h=-0.02, lty=2)
abline(h=0.02, lty=2)
abline(v=0.02, lty=2)
abline(v=0.01, lty=2)

hist(log10(lmis$F_MISS), ylab="Number of SNPs",xlab="Fraction of missing genotypes",main="Fraction of missing data")
#dashed vertical line corresponding to 4% missing data
abline(v=log10(0.05),lty=2)


#close .pdf file
dev.off()



