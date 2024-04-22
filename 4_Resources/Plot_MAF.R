setwd('//lustre/home/sww208/GoDMC/QC/OrganizedSNParray/2_ProcessedData')
#read in frequencies file
freq=read.table("QCoutput/variant_freq.frq", header=T)
#open .png file
pdf("QCoutput/MAF.pdf")
#set output format
par(mfrow=c(1,1))
#plot histogram of frequencies
hist(freq$MAF, ylab="Number of SNPs",xlab="MAF",main="Minor allele frequencies")
#dashed vertical line corresponding to 1% MAF
abline(v=0.02,lty=2)
#close .pdf file
dev.off()