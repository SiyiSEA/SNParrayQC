# plot ibd

arguments <- commandArgs(T)
plinkidb <- arguments[1]
output <- arguments[2]
population <- arguments[3]

#read in caseconpruned.genome file
genome=read.table(plinkidb, header=T, as.is=T)

#The genome file is very big as it contains all possible pairs of
# we look only at the pairs of individuals with PI_HAT > 0.1875, which corresponds to a half way between second # and third degree relatives.

if (dim(genome[genome$PI_HAT > 0.1875,])[1] == 0){
  genome = genome
  message("There is no individual related between second and third degree")
}else{
  genome=genome[genome$PI_HAT > 0.1875,]
} 

#compute Mean(IBD)
mean.ibd=0*genome$Z0 + 1*genome$Z1 + 2*genome$Z2
#compute Var(IBD)
var.ibd=((0 -mean.ibd)^2)*genome$Z0 +
  ((1 -mean.ibd)^2)*genome$Z1 +
  ((2 -mean.ibd)^2)*genome$Z2
#compute SE(IBD)
se.ibd=sqrt(var.ibd)
#open .png file
png(paste0(output, '/', population, "_ibd.png"))
#plot SE(IBD) vs Mean(IBD)
plot(mean.ibd, se.ibd, xlab="Mean IBD", ylab="SE IBD")
#close .png file
dev.off()


duplicate=genome[mean.ibd == 2,]

if (dim(duplicate)[1] == 0 ){
  message("There is no monozygotic tiwns or duplicates.")
}else{
  message(paste("There are ", dim(duplicate)[1], " individuals are monozygotic tiwns either duplicates."))
  #Save the ID of one in each pair into a file fail_ibd_qc.txt for subsequent removal:
  fail_ibd_qc=data.frame(FID=duplicate$FID2, IID=duplicate$IID2)
  write.table(fail_ibd_qc, paste0(output, '/', population, "fail_ibd_qc.txt"), row.names=F, col.names=T, quote=F)
}



