
args<-commandArgs(TRUE)

data_path <- args[1]
Michigan_filename<-args[2]
Sanger_filename<-args[3]
panel <- args[4]

library(data.table)
library(ggplot2)
setwd(data_path)

message("loading the data")
Michiganinfo = fread(Michigan_filename, data.table = F, header = T)
Michiganinfo$Group = paste0("Michigan", panel)

Sangerinfo = fread(Sanger_filename, data.table = F, header = T)
Sangerinfo$Group = paste0("Sanger", panel)


message("plotting the scatter plots and histogram plots")

pdf(paste0(panel,"DesityPlot.pdf"))
par(mfrow=c(2,2))

plot(Michiganinfo$MAF[1:5000], Michiganinfo$R2[1:5000], xlab = "MAF", ylab = "R^2", main = paste0("Michigan", panel))
hist(Michiganinfo$R2, xlab = "R2", main = paste0("Michigan", panel))

plot(Sangerinfo$MAF[1:10000], Sangerinfo$INFO[1:10000], xlab = "MAF", ylab = "Info", main = paste0("Sanger", panel))
hist(Sangerinfo$INFO, xlab = "INFO", main = paste0("Sanger", panel))

par(mfrow=c(1,1))
MAFHRC = rbind(Michiganinfo[,c("MAF", "Group")], Sangerinfo[,c("MAF", "Group")])
ggplot(MAFHRC, aes(x = MAF, fill = Group)) + geom_density(alpha = 0.3) + 
  labs(title = "Kernel Density Plot of MAF between MichiganHRC and SangerHRC",
       x = "MAF",
       y = "Density")
dev.off()

message("The follwing pdf has been generated successfully.")
message(paste0(panel,"DesityPlot.pdf"))