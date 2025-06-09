#PCA total, abiotic - Legionella project.
#All variables numeric

#PCA water total data----
metadata <- read.table("/home/oded/Documents/orit_sivan/articles/oded/PCA/25.04.30/season/sample-metadata.tsv", sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="")
str(metadata)
metadata$moisture_gravimetric <-as.numeric(metadata$moisture_gravimetric)
metadata$depth <-as.factor(metadata$depth)
str(metadata)

#log transformation (+1) - so no inf values+ combine the columns together to one file 
#(the numeric columns in the log file and the factor columns from the original file)
  log.metadata <- log(1+metadata[c(6:14)])
  
  
  
  cbind.log.metadata <-cbind(log.metadata, metadata[c(1:4)])
  str(cbind.log.metadata)
  # now we remove missing values and inf: 
  cbind.log.metadata.NA <- cbind.log.metadata[complete.cases(cbind.log.metadata), ]
  str(cbind.log.metadata.NA)
  

# First we will omit variables with high correlation (>0.7)
  
  #calculate correlation matrix with p value - maybe only for relevant ones - at specific depth/year/stratification
  # bsaed on normalization and n, choose pearson / spearman
  library(Hmisc)
  str(log.metadata)
  cor_rem <- log.metadata[, c(1:9)]
  print(cor_rem)
  cor_pca_Na_conduct <- rcorr(as.matrix(cor_rem), type=c("spearman"))
  print(cor_pca_Na_conduct)
  
  # make the presentation "flattened"
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  cor_total_PCA <- flattenCorrMatrix(cor_pca_Na_conduct$r, cor_pca_Na_conduct$P)
  cor_total_PCA
  
  # Order cor_total_PCA by p-values in ascending order
  cor_total_PCA <- cor_total_PCA[order(cor_total_PCA$p), ]
  
  write.csv(cor_total_PCA,"/home/oded/Documents/orit_sivan/articles/oded/PCA/25.04.30/season/25.04.30_cor_pcoa_total.csv", row.names = FALSE)
  
  #plot the correlation between variables (include histogram)
  library("PerformanceAnalytics")
  library(GGally)
  
    cor_pca_Na_conduct_plot <- cbind.log.metadata.NA[, c(1:9)]
  PcoA.water.correlation <- ggpairs(cor_pca_Na_conduct_plot[1:9], lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.1)))
  PcoA.water.correlation
  
str(metadata)
str(log.metadata)
# normality tests

shapiro <- lapply(log.metadata, shapiro.test)

# almost all variables don't distribute normally.
# this is technically OK for PCR, but skewed data may influence the PCA, so I will center-scale the data in the PCA



# first step is to identify the relevance of all the vectors
# 1. PCA initial vector analysis of all data.

# using prcomp on all the variables, we get the explained variencen for each PC
cbind.log.metadata.NA.pca <- prcomp(cbind.log.metadata.NA[,c(2:9)], center = TRUE,scale = TRUE)
summary(cbind.log.metadata.NA.pca)



# initial graphs to see that the data spreads evenly throught
# the second graph is to look at the loadings. this will tell us what variables influence the PCA the most 
library(ggbiplot)
library(viridis)
#adjust plot margins (otherwise we get en error - Error in plot.new() : figure margins too large)
par(mar = c(1, 1, 1, 1))
plot(cbind.log.metadata.NA.pca$x[,1], cbind.log.metadata.NA.pca$x[,2])
loading.pcoa <- biplot(cbind.log.metadata.NA.pca, scale=0) + theme_classic()
loading.pcoa

pdf("/home/oded/Documents/orit_sivan/articles/oded/PCA/25.04.30/season/25.04.30_season_PcoA_loadings_alaska_chem_all.pdf",  width = 11, height = 8, family="ArialMT")
par(mar = c(1, 1, 1, 1))
loading.pcoa <- biplot(cbind.log.metadata.NA.pca, scale=0) + theme_classic()
dev.off()

loading_pcoa_total <- cbind.log.metadata.NA.pca$rotation[,1:2]
write.csv(loading_pcoa_total,"/home/oded/Documents/orit_sivan/articles/oded/PCA/25.04.30/season/25.04.30_season_loadings_pcoa_total_alaska.csv", row.names = TRUE)

# Now we look at the Rotation or the eigenvectors - these specify the orientation of the principal components relative to the original variables. 
#The elements of an eigenvector, that is, the values within a particular row of matrix A, are the weights aij. These values are called the loadings, 
#and they describe how much each variable contributes to a particular principal component.
#Large loadings (positive or negative) indicate that a particular variable has a strong relationship to a particular principal component
cbind.log.metadata.NA.pca

# in addition we can see which principal componants are relevant. we square the sdev of each
# values > 1 are considered 
PC_sdev <- cbind.log.metadata.NA.pca$sdev^2
PC_sdev

# plot all data
PCA_total <-cbind(cbind.log.metadata.NA, cbind.log.metadata.NA.pca$x[,1:3])
str(PCA_total)
summary(cbind.log.metadata.NA.pca)
cbind.log.metadata.NA.pca

library(ggplot2)
library(viridis)
theme<-theme(panel.background = element_blank(), axis.line = element_line(linewidth = 0.5, colour = "black"), plot.margin=unit(c(1,1,1,1),"line"), legend.position="top", legend.box = "horizontal", legend.key = element_blank())
PcoA <- ggplot(PCA_total, aes(PC1,PC2, colour=exp, shape=core_old, loadings = TRUE))+
  xlab("PC 1 (xxx%)") + ylab("PC 2 (xxx%)") + 
 #geom_point(aes(size = 7))+
geom_label(aes(label= depth, size = 7))+
  scale_colour_manual(values = c(WME = "#3c93c2",SME = "#a6ad4b"))+
  #scale_shape_manual(values = c(9, 16, 8))+
  geom_hline(yintercept=0,linetype="dashed", color = "red", linewidth=0.5, alpha = 0.5) +
  geom_vline(xintercept=0,linetype="dashed", color = "red", linewidth=0.5, alpha = 0.5) +
theme +
  theme(legend.key.size = unit(3,"line")) + labs(colour = "Core")

PcoA

pdf("/home/oded/Documents/orit_sivan/articles/oded/PCA/25.04.30/season/25.04.30_PcoA_season_alaska_chem_all.pdf",  width = 15, height = 8, family="ArialMT")
PcoA
dev.off()

# final graph omitting variables based on loading scores
# using prcomp on all the variables, we get the explained variencen for each PC
cbind.log.metadata.NA.pca <- prcomp(cbind.log.metadata.NA[,c(2:7)], center = TRUE,scale = TRUE)
summary(cbind.log.metadata.NA.pca)

# the second graph is to look at the loading. this will tell us what variables influence the PCA the most 
library(ggbiplot)
library(viridis)
#adjust plot margins (otherwise we get en error - Error in plot.new() : figure margins too large)
par(mar = c(1, 1, 1, 1))
plot(cbind.log.metadata.NA.pca$x[,1], cbind.log.metadata.NA.pca$x[,2])
loading.pcoa <- biplot(cbind.log.metadata.NA.pca, scale=0) + theme_classic()


pdf("/home/oded/Documents/orit_sivan/articles/oded/PCA/25.04.30/season/red/25.04.30_PcoA_red_loadings_chem_all.pdf",  width = 11, height = 8, family="ArialMT")
par(mar = c(1, 1, 1, 1))
loading.pcoa <- biplot(cbind.log.metadata.NA.pca, scale=0) + theme_classic()
dev.off()

loading_pcoa_total <- cbind.log.metadata.NA.pca$rotation[,1:2]
write.csv(loading_pcoa_total,"/home/oded/Documents/orit_sivan/articles/oded/PCA/25.04.30/season/red/25.04.30_red_loadings_pcoa_total.csv", row.names = TRUE)

# Now we look at the Rotation or the eigenvectors - these specify the orientation of the principal components relative to the original variables. 
#The elements of an eigenvector, that is, the values within a particular row of matrix A, are the weights aij. These values are called the loadings, 
#and they describe how much each variable contributes to a particular principal component.
#Large loadings (positive or negative) indicate that a particular variable has a strong relationship to a particular principal component
cbind.log.metadata.NA.pca

# in addition we can see which principal componants are relevant. we square the sdev of each
# values > 1 are considered 
PC_sdev <- cbind.log.metadata.NA.pca$sdev^2
PC_sdev

# plot all data
PCA_total <-cbind(cbind.log.metadata.NA, cbind.log.metadata.NA.pca$x[,1:3])
str(PCA_total)
summary(cbind.log.metadata.NA.pca)
cbind.log.metadata.NA.pca

theme<-theme(panel.background = element_blank(), axis.line = element_line(linewidth = 0.5, colour = "black"), plot.margin=unit(c(1,1,1,1),"line"), legend.position="top", legend.box = "horizontal", legend.key = element_blank())
PcoA_red <- ggplot(PCA_total, aes(PC1,PC2, colour=exp, shape=core_old, loadings = TRUE))+
  xlab("PC 1 (37.82%)") + ylab("PC 2 (25.44%)") + 
  #geom_point(aes(size = 7))+
  geom_label(aes(label=depth, size = 7))+
  scale_colour_manual(values = c(WME = "#3c93c2",SME = "#a6ad4b"))+
  #scale_shape_manual(values = c(9, 16, 8))+
  geom_hline(yintercept=0,linetype="dashed", color = "red", linewidth=0.5, alpha = 0.5) +
  geom_vline(xintercept=0,linetype="dashed", color = "red", linewidth=0.5, alpha = 0.5) +
  theme +
  theme(legend.key.size = unit(3,"line")) + labs(colour = "Layer", shape = "Mixis")

PcoA_red

pdf("/home/oded/Documents/orit_sivan/articles/oded/PCA/25.04.30/season/red/25.04.30_red_PcoA_chem_all.pdf",  width = 15, height = 8, family="ArialMT")
PcoA_red
dev.off()
