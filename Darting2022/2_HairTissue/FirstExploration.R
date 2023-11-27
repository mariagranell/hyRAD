

#################################################

setwd("/Users/mariagranell/Repositories/hyRAD/Darting2022/2_HairTissue")

#################################################

#install.packages("HardyWeinberg")
library(mice)
library(HardyWeinberg)
library(gaston)
library(SNPRelate)


add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


############### READ THE DATA #################

# Reads the VCF
#VCF_Vervet = gaston::read.vcf("populations.snps.vcf", convert.chr = FALSE, get.info = TRUE) #

# Transfrom VCF to gds (run only once)
#DATA = snpgdsVCF2GDS("populations.snps.vcf", "populations.snps.gds")
# and then open the gds
#DATA = snpgdsOpen("populations.snps.gds")

# to close the GDS run :
#showfile.gds(closeall=TRUE)

# read the gds as a doseage matrix
gen_mat <- snpgdsGetGeno('populations.snps.gds')

############### EXPLORE THE DATSET ##############

# Visualise the misigness per site
hist(colSums(is.na(gen_mat))/93)

# Visualise the misigness per indiv
hist(rowSums(is.na(gen_mat))/ncol(gen_mat), breaks = 100)
# identify individuals with more than 90% missing data
VCF_Vervet@ped$id[rowSums(is.na(gen_mat))/ncol(gen_mat)>0.9]

# Nombre d'individus heterozygote : colSums(gen_mat==1, na.rm = TRUE)
# nombre d'individus genotypes : colSums(! is.na(gen_mat))
Ho = colSums(gen_mat==1, na.rm = TRUE)/colSums(! is.na(gen_mat))

# 
plot(VCF_Vervet@p, Ho, pch = 20, col = add.alpha("blue", 0.2))
lines(seq(0,0.5,0.01), 2*seq(0,0.5,0.01)*(1-seq(0,0.5,0.01)))

####

# reduc the vcf for SNPs gentoyped in more than 60 ind, because less that that is not really usefull
VCF_Vervet_Reduc = VCF_Vervet[,colSums(! is.na(gen_mat))>=60]

hw.mpVervet <- HardyWeinberg::HWExactStats(cbind(VCF_Vervet_Reduc@snps$N0, VCF_Vervet_Reduc@snps$N1, VCF_Vervet_Reduc@snps$N2), midp=TRUE)

outliers <- which(-log10(hw.mpVervet) > 6)
plot(VCF_Vervet_Reduc@p, VCF_Vervet_Reduc@snps$hz, col="black", pch=16, cex=0.6, xlab='p', ylab='Het')
points(VCF_Vervet_Reduc[,outliers]@p, VCF_Vervet_Reduc[,outliers]@snps$hz, col="red",pch=16,cex=0.6)
lines(seq(0,0.5,0.01), 2*seq(0,0.5,0.01)*(1-seq(0,0.5,0.01)), col = "blue")

##############################################

TissueType = substr(VCF_Vervet@ped$id, 1, 2)
ID = substr(VCF_Vervet@ped$id, 4, 5)

PCA_All = snpgdsPCA(DATA, autosome.only = FALSE, verbose = TRUE, eigen.cnt = 0, remove.monosnp=TRUE, maf = 0.05, missing.rate = 0.2)
plot(PCA_All$eigenvect[,1],PCA_All$eigenvect[,2], pch = 20, col = c("grey20", "grey80")[as.numeric(as.factor(TissueType))])
text(PCA_All$eigenvect[,1],PCA_All$eigenvect[,2], ID, cex = 0.6)

library(ade4)
plot(PCA_All$eigenvect[,1],PCA_All$eigenvect[,2], pch = 20, col = c("grey20", "grey80")[as.numeric(as.factor(TissueType))])
s.class(PCA_All$eigenvect[,1:2], fac = as.factor(ID), grid = FALSE, clabel = 0, cpoint = 0, cellipse = 0, add.plot = TRUE)


################################################

SNPsID = snpgdsSNPList(DATA)
SNPsID[GoodSNPs]

ImbBeta = snpgdsIndivBeta(DATA, snp.id = SNPsID$snp.id[GoodSNPs], sample.id = VCF_Vervet@ped$id[TissueType=="DT"] ,autosome.only = FALSE, inbreeding = TRUE)
heatmap(ImbBeta$beta, scale = "none", cexRow = 0.4, cexCol = 0.4)

ImbBeta = snpgdsIndivBeta(DATA, snp.id = SNPsID$snp.id[GoodSNPs], sample.id = VCF_Vervet@ped$id[TissueType=="DT"] ,autosome.only = FALSE, inbreeding = TRUE)
heatmap(ImbBeta$beta, scale = "none", cexRow = 0.4, cexCol = 0.4)

hist(diag(ImbBeta$beta))
hist(ImbBeta$beta[upper.tri(ImbBeta$beta)])

ImbBeta = snpgdsIndivBeta(DATA, snp.id = SNPsID$snp.id[GoodSNPs], sample.id = VCF_Vervet@ped$id[TissueType=="DT" & !(ID %in% c(3,13,23))] ,autosome.only = FALSE, inbreeding = TRUE)
heatmap(ImbBeta$beta, scale = "none", cexRow = 0.4, cexCol = 0.4)
BetaMat = ImbBeta$beta
diag(BetaMat) = NA
colnames(BetaMat) = VCF_Vervet@ped$id[TissueType=="DT" & !(ID %in% c("3","13","23", "11", "9-"))]
rownames(BetaMat) = VCF_Vervet@ped$id[TissueType=="DT" & !(ID %in% c(3,13,23))]
heatmap(BetaMat, scale = "none", cexRow = 0.4, cexCol = 0.4)

hist(diag(ImbBeta$beta), breaks = 30)
hist(ImbBeta$beta[upper.tri(ImbBeta$beta)])

hist(BetaMat[upper.tri(BetaMat)], breaks = 100)

#########################################################################

ImbBeta = snpgdsIndivBeta(DATA, snp.id = SNPsID$snp.id[GoodSNPs], sample.id = VCF_Vervet@ped$id[TissueType=="DT" & !(ID %in% c("3","13","23", "11", "9-", "20"))] ,autosome.only = FALSE, inbreeding = TRUE)
heatmap(ImbBeta$beta, scale = "none", cexRow = 0.4, cexCol = 0.4)
BetaMat = ImbBeta$beta
diag(BetaMat) = NA
colnames(BetaMat) = VCF_Vervet@ped$id[TissueType=="DT" & !(ID %in% c("3","13","23", "11", "9-", "20"))]
rownames(BetaMat) = VCF_Vervet@ped$id[TissueType=="DT" & !(ID %in% c("3","13","23", "11", "9-", "20"))]
heatmap(BetaMat, scale = "none", cexRow = 0.4, cexCol = 0.4)

hist(BetaMat[upper.tri(BetaMat)], breaks = 100)
