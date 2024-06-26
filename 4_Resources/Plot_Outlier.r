##########################################
# Color all the outliers for each dataset based on PC1 and PC2
library(ggplot2)
ColorOutliers <- function(eigenFileName, outlierFileName) { 
  
  outname<-unlist(strsplit(eigenFileName, "[.]"))[[1]]
  prefix <- unlist(strsplit(outname, "[/]"))[[9]]
  
  data = read.table(eigenFileName)
  outlier = read.table(outlierFileName, header = T)
  
  data$Group = "Sample"
  OutlierPC1 = outlier[which(outlier$nPC == "PC1"),"FID"]
  OutlierPC2 = outlier[which(outlier$nPC == "PC2"),"FID"]
  data[which(data$V1 %in% OutlierPC1) ,'Group'] = "PC1Outliers"
  data[which(data$V1 %in% OutlierPC2) ,'Group'] = "PC2Outliers"
  
  ggplot(data, aes(x = V3, y = V4, color = Group)) +
    geom_point() + 
    xlab("PC1") + 
    ylab("PC2") +
    labs(title = paste(prefix,"with",length(unique(c(OutlierPC1, OutlierPC2))), "Outliers")) 
}

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/scz_ab_eur_QCd.imqc.pca.eigenvec"
outlierFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/scz_ab_eur_QCd_OutliersFromPC_3SDfromMean_temp.txt"
ColorOutliers(eigenFileName, outlierFileName)

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Michigan.imqc.pca.eigenvec"
outlierFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Michigan_OutliersFromPC_3SDfromMean_temp.txt"
ColorOutliers(eigenFileName, outlierFileName)

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Michigan.imqc.pca.eigenvec"
outlierFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Michigan_OutliersFromPC_3SDfromMean_temp.txt"
ColorOutliers(eigenFileName, outlierFileName)

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Sanger.imqc.pca.eigenvec"
outlierFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Sanger_OutliersFromPC_3SDfromMean_temp.txt"
ColorOutliers(eigenFileName, outlierFileName)

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Sanger.imqc.pca.eigenvec"
outlierFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Sanger_OutliersFromPC_3SDfromMean_temp.txt"
ColorOutliers(eigenFileName, outlierFileName)

##############################################
# QCd outliers in Sanger and Michigan
ColorOutliers <- function(eigenFileName, outlierFileName) { 
  
  outname<-unlist(strsplit(eigenFileName, "[.]"))[[1]]
  prefix <- unlist(strsplit(outname, "[/]"))[[9]]
  
  data = read.table(eigenFileName)
  outlier = read.table(outlierFileName, header = T)
  
  data$Group = "Sample"
  OutlierPC1 = outlier[which(outlier$nPC == "PC1"),"FID"]
  OutlierPC2 = outlier[which(outlier$nPC == "PC2"),"FID"]
  data[which(data$V1 %in% OutlierPC1) ,'Group'] = "PC1Outliers"
  data[which(data$V1 %in% OutlierPC2) ,'Group'] = "PC2Outliers"
  
  ggplot(data, aes(x = V3, y = V4, color = Group)) +
    geom_point() + 
    xlab("PC1") + 
    ylab("PC2") +
    labs(title = paste(prefix,"with QCd",length(unique(c(OutlierPC1, OutlierPC2))), "Outliers")) 
}

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Michigan.imqc.pca.eigenvec"
outlierFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/scz_ab_eur_QCd_OutliersFromPC_3SDfromMean_temp.txt"
ColorOutliers(eigenFileName, outlierFileName)

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Michigan.imqc.pca.eigenvec"
outlierFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/scz_ab_eur_QCd_OutliersFromPC_3SDfromMean_temp.txt"
ColorOutliers(eigenFileName, outlierFileName)

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Sanger.imqc.pca.eigenvec"
outlierFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/scz_ab_eur_QCd_OutliersFromPC_3SDfromMean_temp.txt"
ColorOutliers(eigenFileName, outlierFileName)

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Sanger.imqc.pca.eigenvec"
outlierFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/scz_ab_eur_QCd_OutliersFromPC_3SDfromMean_temp.txt"
ColorOutliers(eigenFileName, outlierFileName)

##############################################
# Sanger outliers in QCd
ColorOutliers <- function(eigenFileName, outlierFileName, outliername) { 
  
  outname<-unlist(strsplit(eigenFileName, "[.]"))[[1]]
  prefix <- unlist(strsplit(outname, "[/]"))[[9]]
  
  data = read.table(eigenFileName)
  outlier = read.table(outlierFileName, header = T)
  
  data$Group = "Sample"
  OutlierPC1 = outlier[which(outlier$nPC == "PC1"),"FID"]
  OutlierPC2 = outlier[which(outlier$nPC == "PC2"),"FID"]
  data[which(data$V1 %in% OutlierPC1) ,'Group'] = "PC1Outliers"
  data[which(data$V1 %in% OutlierPC2) ,'Group'] = "PC2Outliers"
  
  ggplot(data, aes(x = V3, y = V4, color = Group)) +
    geom_point() + 
    xlab("PC1") + 
    ylab("PC2") +
    labs(title = paste(prefix,"with",outliername,length(unique(c(OutlierPC1, OutlierPC2))), "Outliers")) 
}

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/scz_ab_eur_QCd.imqc.pca.eigenvec"
outlierFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Michigan_OutliersFromPC_3SDfromMean_temp.txt"
ColorOutliers(eigenFileName, outlierFileName, outliername = "Michigan1000G")

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/scz_ab_eur_QCd.imqc.pca.eigenvec"
outlierFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Michigan_OutliersFromPC_3SDfromMean_temp.txt"
ColorOutliers(eigenFileName, outlierFileName, "MichiganHRC")


eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/scz_ab_eur_QCd.imqc.pca.eigenvec"
outlierFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Sanger_OutliersFromPC_3SDfromMean_temp.txt"
ColorOutliers(eigenFileName, outlierFileName, "SangerHRC")

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/scz_ab_eur_QCd.imqc.pca.eigenvec"
outlierFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Sanger_OutliersFromPC_3SDfromMean_temp.txt"
ColorOutliers(eigenFileName, outlierFileName, "Sanger1000G")

##############################################
# Michigan outliers in QCd

MichiganOutlier = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Michigan_OutliersFromPC_3SDfromMean_temp.txt"
MichiganOutlier = read.table(MichiganOutlier, header = T)

SangerOutlier = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Sanger_OutliersFromPC_3SDfromMean_temp.txt"
SangerOutlier = read.table(SangerOutlier, header = T)

intersect(MichiganOutlier[which(MichiganOutlier$nPC == "PC1"),'FID'], SangerOutlier[which(SangerOutlier$nPC == "PC1"),'FID'])
length(intersect(MichiganOutlier[which(MichiganOutlier$nPC == "PC1"),'FID'], SangerOutlier[which(SangerOutlier$nPC == "PC1"),'FID'])
)

length(intersect(MichiganOutlier[which(MichiganOutlier$nPC == "PC2"),'FID'], SangerOutlier[which(SangerOutlier$nPC == "PC2"),'FID'])
)

MichiganOutlier = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Michigan_OutliersFromPC_3SDfromMean_temp.txt"
MichiganOutlier = read.table(MichiganOutlier, header = T)

SangerOutlier = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Sanger_OutliersFromPC_3SDfromMean_temp.txt"
SangerOutlier = read.table(SangerOutlier, header = T)

intersect(MichiganOutlier[which(MichiganOutlier$nPC == "PC1"),'FID'], SangerOutlier[which(SangerOutlier$nPC == "PC1"),'FID'])
length(intersect(MichiganOutlier[which(MichiganOutlier$nPC == "PC1"),'FID'], SangerOutlier[which(SangerOutlier$nPC == "PC1"),'FID'])
)

length(intersect(MichiganOutlier[which(MichiganOutlier$nPC == "PC2"),'FID'], SangerOutlier[which(SangerOutlier$nPC == "PC2"),'FID'])
)

####################################################
# plot
eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Sanger.imqc.pca.eigenvec"
outlierFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Michigan_OutliersFromPC_3SDfromMean_temp.txt"
ColorOutliers(eigenFileName, outlierFileName, "MichiganHRC")

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Michigan.imqc.pca.eigenvec"
outlierFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Sanger_OutliersFromPC_3SDfromMean_temp.txt"
ColorOutliers(eigenFileName, outlierFileName, "SangerHRC")

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Sanger.imqc.pca.eigenvec"
outlierFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Michigan_OutliersFromPC_3SDfromMean_temp.txt"
ColorOutliers(eigenFileName, outlierFileName, "Michigan1000G")

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Michigan.imqc.pca.eigenvec"
outlierFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Sanger_OutliersFromPC_3SDfromMean_temp.txt"
ColorOutliers(eigenFileName, outlierFileName, "Sanger1000G")

####################################################
# SD

CalcualteSD <- function(eigenFileName) { 
  
  outname<-unlist(strsplit(eigenFileName, "[.]"))[[1]]
  prefix <- unlist(strsplit(outname, "[/]"))[[9]]
  
  data = read.table(eigenFileName)
  PC1sd = sd(data$V3)
  PC2sd = sd(data$V4)
  print(paste(prefix, "has SD of", PC1var, "for PC1."))
  print(paste(prefix, "has SD of", PC2var, "for PC2."))
  print(paste(prefix, "has variance of", PC1var^2, "for PC1."))
  print(paste(prefix, "has variance of", PC2var^2, "for PC2."))
}

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Michigan.imqc.pca.eigenvec"
CalcualteSD(eigenFileName)

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Michigan.imqc.pca.eigenvec"
CalcualteSD(eigenFileName)

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Sanger.imqc.pca.eigenvec"
CalcualteSD(eigenFileName)

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Sanger.imqc.pca.eigenvec"
CalcualteSD(eigenFileName)

eigenFileName = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/scz_ab_eur_QCd.imqc.pca.eigenvec"
CalcualteSD(eigenFileName)

####################################################
# PC1 vs PC1
RelationPCs <- function(eigenFileName1, eigenFileName2){
  outname1<-unlist(strsplit(eigenFileName1, "[.]"))[[1]]
  prefix1 <- unlist(strsplit(outname1, "[/]"))[[9]]
  
  outname2<-unlist(strsplit(eigenFileName2, "[.]"))[[1]]
  prefix2 <- unlist(strsplit(outname2, "[/]"))[[9]]
  
  data1 = read.table(eigenFileName1)
  data2 = read.table(eigenFileName2)
  
  plot(data1$V3, data1$V3, xlab = paste(prefix1,"PC1"), ylab=paste(prefix2, "PC1"))
}

eigenFileName1 = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Sanger.imqc.pca.eigenvec"
eigenFileName2 = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Michigan.imqc.pca.eigenvec"
RelationPCs(eigenFileName1, eigenFileName2)

eigenFileName1 = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Sanger.imqc.pca.eigenvec"
eigenFileName2 = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Michigan.imqc.pca.eigenvec"
RelationPCs(eigenFileName1, eigenFileName2)

eigenFileName1 = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/scz_ab_eur_QCd.imqc.pca.eigenvec"
eigenFileName2 = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Michigan.imqc.pca.eigenvec"
RelationPCs(eigenFileName1, eigenFileName2)

eigenFileName1 = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/scz_ab_eur_QCd.imqc.pca.eigenvec"
eigenFileName2 = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_HRC_Sanger.imqc.pca.eigenvec"
RelationPCs(eigenFileName1, eigenFileName2)

eigenFileName1 = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/scz_ab_eur_QCd.imqc.pca.eigenvec"
eigenFileName2 = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Sanger.imqc.pca.eigenvec"
RelationPCs(eigenFileName1, eigenFileName2)

eigenFileName1 = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/scz_ab_eur_QCd.imqc.pca.eigenvec"
eigenFileName2 = "/lustre/home/sww208/QC/SNParrayQC/3_Results/PCAVariants/data_filtered_1000G_Michigan.imqc.pca.eigenvec"
RelationPCs(eigenFileName1, eigenFileName2)
