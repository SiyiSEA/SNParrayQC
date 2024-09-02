suppressPackageStartupMessages(library(bigsnpr))
suppressPackageStartupMessages(library(ggplot2))

arguments <- commandArgs(T)
pathtodir <- arguments[1]
bedfile <- arguments[2]
dataname <- arguments[3]
Sthreshold <- as.numeric(arguments[4])
message("This is S threshold ", Sthreshold)
homothreshold <- as.numeric(arguments[5])
message("This is homothreshold ", homothreshold)

setwd(pathtodir)

# read in the PLINK files
snpfiles = snp_readBed(bedfile, backingfile = sub_bed(bedfile))
# Attach the "bigSNP" object in R session
snpattach = snp_attach(snpfiles)

# see how it looks like
obj.bed <- bed(bedfile)
obj.bed

# first time to detect the long-range LD with inital clumping r^2>0.2
# in this case, the long-range LD won't affect the PCA
obj.svd <- bed_autoSVD(obj.bed, k = 20, ncores = nb_cores())


# detect the outlier samples in PCA
prob <- bigutilsr::prob_dist(obj.svd$u, ncores = nb_cores())
S <- prob$dist.self / sqrt(prob$dist.nn)

# histgram for select the threshold
png(paste0("hist_SScore", dataname, ".png"))
hist(S, breaks = "FD", xlab = "Statistic of outlierness",
     ylab = "Frequency (sqrt-scale)",
     main = "Distribution of statistics (S)")
dev.off()

# PCA for select the threshold
PCA_SScore = plot_grid(plotlist = lapply(1:10, function(k) {
  plot(obj.svd, type = "scores", scores = 2 * k - 1:0, coeff = 0.6) +
    aes(color = S) +
    scale_colour_viridis_c()
}), scale = 0.95)
ggsave(paste0("PCA_SScore", dataname, ".pdf"), plot=PCA_SScore, width = 16, height = 15)


if (is.na(Sthreshold)){
  message("Please detect the threshold based on the hist and PCA plot with S score.")
  message("Skip the removing outliers by S score.")
  PC <- predict(obj.svd)
  ind.norel <- rows_along(obj.bed)
  ind.row <- ind.norel
  }else{
  PCA_outliers_SScore = plot_grid(plotlist = lapply(1:10, function(k) {
    plot(obj.svd, type = "scores", scores = 2 * k - 1:0, coeff = 0.6) +
      aes(color = S > as.numeric(Sthreshold)) +  # threshold based on histogram
      scale_colour_viridis_d()
  }), scale = 0.95)
  ggsave(paste0("PCA_outliers_SScore", dataname, ".pdf"), plot=PCA_outliers_SScore, width = 16, height = 15)

  # remove the outliers
  message("Removing outliers identified by the S score.")
  ind.norel = rows_along(obj.bed)
  ind.row <- ind.norel[S < as.numeric(Sthreshold)]
  ind.col <- attr(obj.svd, "subset")
  obj.svd2 <- bed_autoSVD(obj.bed, ind.row = ind.row,
                          ind.col = ind.col, thr.r2 = NA,
                          k = 20, ncores = nb_cores())

  # plot(obj.svd2, type = "loadings", loadings = 1:20, coeff = 0.4)
  remove_outliers_SScore = plot(obj.svd2, type = "scores", scores = 1:20, coeff = 0.4)
  ggsave(paste0("PCA_remove_outliers_SScore", dataname, ".pdf"), plot=remove_outliers_SScore, width = 16, height = 15)

  PC <- predict(obj.svd2)
}

# dectect samples have a different ancestry from most of the samples in the data
ldist <- log(bigutilsr::dist_ogk(PC))

# histgram for select the threshold
png(paste0("hist_Mandist", dataname, ".png"))
hist(ldist, breaks = "FD", xlab = "log Mahalanobis distance", main = "Distribution of (log squared) distances")
dev.off()

# PCA for select the threshold
PCA_Mandist = plot_grid(plotlist = lapply(1:9, function(k) {
  k1 <- 2 * k - 1; k2 <- 2 * k
  qplot(PC[, k1], PC[, k2]) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2)) +
    theme_bigstatsr(0.5) +
    aes(color = ldist) +
    theme(legend.position = "left") +
    scale_color_viridis_c(trans = "log") +
    coord_equal()
}), ncol = 3)
ggsave(paste0("PCA_Mandist", dataname, ".pdf"), plot=PCA_Mandist, width = 16, height = 15)

if (is.na(homothreshold)){
  message("Please detect the threshold based on hist and PCA plot with Mahalanobis distance in log scale!")
  message("Skip the removing outliers identified by the homogeneous")
  homo.row <- ind.row
  }else{
  PCA_outliers_Mandist = plot_grid(plotlist = lapply(1:9, function(k) {
    k1 <- 2 * k - 1; k2 <- 2 * k
    qplot(PC[, k1], PC[, k2]) +
      labs(x = paste0("PC", k1), y = paste0("PC", k2)) +
      theme_bigstatsr(0.7) +
      aes(color = ldist < as.numeric(homothreshold)) +
      scale_color_viridis_d(direction = -1) +
      coord_equal()
  }), ncol = 3)
  ggsave(paste0("PCA_outliers_Mandist", dataname, ".pdf"), plot=PCA_outliers_Mandist, width = 16, height = 15)

  # remove the outliers
  message("Removing outliers identified by the homogeneous.")
  homo.row <- ind.row[ldist < homothreshold]
  homo.col <- attr(obj.svd2, "subset")
  obj.svd3 <- bed_autoSVD(obj.bed, ind.row = homo.row,
                          ind.col = ind.col, thr.r2 = NA,
                          k = 20, ncores = nb_cores())
  # plot(obj.svd3, type = "loadings", loadings = 1:10, coeff = 0.4)
  PCA_remove_outliers_Mandist = plot(obj.svd3, type = "scores", scores = 1:10, coeff = 0.4)
  ggsave(paste0("PCA_remove_outliers_Mandist", dataname, ".pdf"), plot=PCA_remove_outliers_Mandist, width = 16, height = 15)

}

# generate a list of sample ID for keeping the remining sample
famfile = paste0(sub_bed(bedfile), ".fam")
fam = read.table(famfile, header = F)
length(homo.row)
keeplist = fam[homo.row, 1:2]
write.table(keeplist, file = paste0(sub_bed(bedfile), ".keep"), sep = " ", quote = F, row.names = F, col.names = F)