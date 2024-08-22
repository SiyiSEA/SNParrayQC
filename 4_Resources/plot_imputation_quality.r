
args<-commandArgs(TRUE)

data_path <- args[1]
Michigan_filename<-args[2]
Sanger_filename<-args[3]
panel <- args[4]


MichiganPanel = T
SangerPanel = T
if (Michigan_filename == "NULL"){
  MichiganPanel = F
}
if (Sanger_filename == "NULL"){
  SangerPanel = F
}

library(data.table)
library(ggplot2)
setwd(data_path)

message("loading the data")

if (MichiganPanel == T){
  message(Michigan_filename)
  Michiganinfo = fread(Michigan_filename, data.table = F, header = T)
  Michiganinfo$Group = paste0("Michigan", panel)
}


if (SangerPanel == T){
  message(Sanger_filename)
  Sangerinfo = fread(Sanger_filename, data.table = F, header = T)
  Sangerinfo$Group = paste0("Sanger", panel)
}


message("plotting the scatter plots and histogram plots")

pdf(paste0(panel,"MAFDesityPlot.pdf"))
par(mfrow=c(1,2))

if (MichiganPanel == T){
  plot(Michiganinfo$MAF[1:5000], Michiganinfo$R2[1:5000], xlab = "MAF", ylab = "Info", main = paste0("Michigan", panel))
  hist(Michiganinfo$R2, xlab = "R2", main = paste0("Michigan", panel))
}

if (SangerPanel == T){
  plot(Sangerinfo$MAF[1:10000], Sangerinfo$INFO[1:10000], xlab = "MAF", ylab = "Info", main = paste0("Sanger", panel))
  hist(Sangerinfo$INFO, xlab = "INFO", main = paste0("Sanger", panel))
}


if (MichiganPanel == T & SangerPanel == T){
  par(mfrow=c(1,1))
  MAFHRC = rbind(Michiganinfo[,c("MAF", "Group")], Sangerinfo[,c("MAF", "Group")])
  ggplot(MAFHRC, aes(x = MAF, fill = Group)) + geom_density(alpha = 0.3) + 
    labs(title = "Kernel Density Plot of MAF between Sanger and Michigan",
       x = "MAF",
       y = "Density")
}



dev.off()

message("The follwing pdf has been generated successfully.")
message(paste0(panel,"MAFDesityPlot.pdf"))
