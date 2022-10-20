

#######################################################################################################
#######################################################################################################
#                   Phenomic Selection Practicals 2022/09/20
# Genomic and Phenomic selection in bread wheat using NIRS measured on grains (Rincent et al. 2018)
# renaud.rincent@inrae.fr and vincent.segura@inrae.fr
#######################################################################################################
#######################################################################################################

rm(list=ls())

# Set your working directory
setwd("D:/sauvegarde/Congres/2022_OrganisationEucarpia/SateliteMeeting_PhenomicSelection/Practical/")

####################################################################################
####################################################################################
# I/ Load packages and data
####################################################################################
####################################################################################

# Load packages
###############################

#install.packages(c("prospectr", "signal", "rrBLUP"))

library(prospectr)
library(signal)
library(rrBLUP)


# Load data
###############################

# Raw NIRS data (measured on grains in the reference environment "EnvRef") :
data_raw <- readRDS("NIRS_Dry.Rds")   # Grain Dry

head(data_raw$NIRS[, 1:5])

# Phenotypic data (Grain Yield) :
pheno <- read.table("Adjmeans_Final.csv", header=T, check.names=FALSE)

head(pheno)

# Genotypic data (SNP) :
#geno <- read.table("GenotypicData.csv", check.names=FALSE)          # full genotyping dataset (84259 markers)
geno <- read.table("GenotypicData_subset.csv", check.names=FALSE)  # Load subset of 10533 markers if your computer is too slow

dim(geno)
geno[1:5, 1:5]

####################################################################################
####################################################################################
# II/ Filters on genomic data and computation of kinship matrix
####################################################################################
####################################################################################

# Filter on MAF
###############################

p <- rowMeans(geno) # average frequency of the reference allele
summary(p)

ToRemove <- which(p <= 0.025 | p >= 0.975) # Remove markers with a MAF below 5%
length(ToRemove)

geno <- geno[-ToRemove, ]
dim(geno)

# Compute Kinship (matA1)
###############################

p <- rowMeans(geno)
q <- 1-p

genot.ok <- 2*t(geno)
rm(geno)

genot.scaled <- scale(genot.ok, center=2*p, scale=sqrt(4*p*q))

matA1 <- tcrossprod(genot.scaled) / ncol(genot.scaled)    # tcrossprod(X) is equivalent to X %*% t(X)
                                              # matA1 is your genomic kinship

rm(p, q, genot.scaled, ToRemove)

####################################################################################
####################################################################################
# III/ Statistical preprocessing of the spectra
####################################################################################
####################################################################################

#Graphical representation of raw spectra:

#data_raw=data.frame(lambda=seq(400,2400),NIRS=I(NIRSmatrix))

matplot(data_raw$lambda, data_raw$NIRS, type="l", lty=1, pch=0,
        xlab = "Lambda (nm)", ylab="Absorbance", xlim=c(390, 2500))

# Normalization of the spectra
###############################

data_norm <- data_raw

data_norm$NIRS <- scale(data_raw$NIRS)

matplot(data_norm$lambda, data_norm$NIRS, type="l", lty=1, pch=0,
        xlab = "Lambda (nm)", ylab="Absorbance", xlim=c(390, 2500))

# Detrend
###############################
# De-trending is performed through subtraction of a linear or polynomial fit of baseline from the original spectrum

data_dt <- data_raw

data_dt$NIRS <- t(detrend(X = t(data_raw$NIRS), wav = data_raw$lambda)) # Standard Normal Variate followed by fitting a 2nd order linear model, 
                                                                            # output are the fitted residuals

matplot(data_dt$lambda, data_dt$NIRS, type="l", lty=1, pch=0,
        xlab = "Lambda (nm)", ylab="Absorbance", xlim=c(390, 2500))


# 1st and 2nd derivatives
###############################

tsf <- (max(data_raw$lambda) - min(data_raw$lambda)) / (length(data_raw$lambda) - 1)
tsf #resolution of spectral data in nm

#1st derivative on raw spectra:

data_der1 <- data_raw

data_der1$NIRS <- as.matrix(apply(data_raw$NIRS, 2, function(x) {
  sgolayfilt(x, p = 2, n = 37, m = 1, ts = tsf)}))   # p is the filter order, m-th derivative, n support points

rownames(data_der1$NIRS) <- rownames(data_raw$NIRS)

matplot(data_der1$lambda, data_der1$NIRS, type="l", lty=1, pch=0,
        xlab = "Lambda (nm)", ylab="Absorbance", xlim=c(390, 2500))


#2nd derivative on raw spectra:

data_der2 <- data_raw
data_der2$NIRS <- as.matrix(apply(data_raw$NIRS, 2, function(x) {
  sgolayfilt(x, p = 3, n = 61, m = 2, ts = tsf)}))

rownames(data_der2$NIRS) <- rownames(data_raw$NIRS)

matplot(data_der2$lambda, data_der2$NIRS, type="l", lty=1, pch=0,
        xlab = "Lambda (nm)", ylab="Absorbance", xlim=c(390, 2500))

#1st derivative on normalized spectra:

data_norm_der1 <- data_norm
data_norm_der1$NIRS <- as.matrix(apply(data_norm$NIRS, 2, function(x) {
  sgolayfilt(x, p = 2, n = 37, m = 1, ts = tsf) }))

rownames(data_norm_der1$NIRS) <- rownames(data_raw$NIRS)

matplot(data_norm_der1$lambda, data_norm_der1$NIRS, type="l", lty=1, pch=0,
        xlab = "Lambda (nm)", ylab="Absorbance", xlim=c(390, 2500))

#2nd derivative on normalized spectra:

data_norm_der2 <- data_norm
data_norm_der2$NIRS <- as.matrix(apply(data_norm$NIRS, 2, function(x) {
  sgolayfilt(x, p = 3, n = 61, m = 2, ts = tsf) }))

rownames(data_norm_der2$NIRS) <- rownames(data_raw$NIRS)

matplot(data_norm_der2$lambda, data_norm_der2$NIRS, type="l", lty=1, pch=0,
        xlab = "Lambda (nm)", ylab="Absorbance", xlim=c(390, 2500))


# Stacking all spectra pre-treatments into a list:
###############################

spectra <- list("raw" = data_raw, "norm" = data_norm, "dt" = data_dt,
                 "der1" = data_der1, "der2" = data_der2,
                 "norm_der1" = data_norm_der1, "norm_der2" = data_norm_der2)

rm(data_raw, data_norm, data_dt, data_der1, data_der2, data_norm_der1, data_norm_der2, tsf)

####################################################################################
####################################################################################
# IV/ Are the spectra under genetic determinism ?
####################################################################################
####################################################################################


# Fit a GBLUP model to each wavelength to estimate the genomic and residual variances
# along the spectrum
###############################

spec <- spectra$norm_der1$NIRS

GenomicVariance <- ResidualVariance <- rep(NA, nrow(spec))

for (i in 1:nrow(spec)) {
  print(i)

  mod4 <- mixed.solve(y = spec[i,], K = matA1)
  GenomicVariance[i] <- mod4$Vu
  ResidualVariance[i] <- mod4$Ve
  rm(mod4)
}

PropGenomicVariance <- GenomicVariance/(GenomicVariance+ResidualVariance)*100

#Graphical representation of the proportion of variance explained by genomics along the spectrum

par(mar=c(4, 4, 4, 4))
plot(seq(400, 2498, by=2), PropGenomicVariance, type="l", xlab="lambda (nm)",
     xlim=c(400, 3000), ylim=c(0, 100), ylab="Proportion of of variance explained by genomics")
polygon(c(400, seq(400, 2498, by=2), 2498), c(0, PropGenomicVariance, 0), col = "brown1")
polygon(c(400, seq(400, 2498, by=2), 2498), c(100, PropGenomicVariance, 100), col = "dodgerblue4")
legend(2500, 90, c("Residual", "Genomic"), lty = 0, bty = "n", fill = c("dodgerblue4", "brown1"), cex=1)

rm(spec, GenomicVariance, ResidualVariance, PropGenomicVariance)

####################################################################################
####################################################################################
# V/ Genomic and Phenomic predictions (within environment cross validations)
####################################################################################
####################################################################################

# Here we use GBLUP but of course you could use any GS model

spct <- spectra$norm_der1$NIRS          # Choose a pretreatment
spct2 <- scale(t(spct), center=T, scale=T) # scale absorbance at each wavelength (predictor)
matH <- tcrossprod(spct2)/ncol(spct2)     # Compute the hyperspectral similarity matrix

Nenvt=8   # Number of environments
Nind=nrow(pheno) # Number of varieties

Nrep=25   # Number of repetition for the cross validation
Nout=30   # Number of varieties in the predicted set


AccuHBLUP <- AccuGBLUP <- matrix(NA, Nrep, Nenvt)
colnames(AccuHBLUP) <- colnames(AccuGBLUP) <- colnames(pheno)[2:ncol(pheno)]


for (envt in 2:9) {
  print(envt)
  phenotype <- pheno[, envt]
  
  for (rep in 1:Nrep) {
    
    valid <- sample(Nind, Nout)
    phenoTrain <- phenotype
    phenoTrain[valid] <- NA
    
    gblup <- mixed.solve(y=phenoTrain, K=matA1)
    hblup <- mixed.solve(y=phenoTrain, K=matH)
    AccuGBLUP[rep,(envt-1)] <- cor(gblup$u[valid], phenotype[valid], use="complete.obs")
    AccuHBLUP[rep,(envt-1)] <- cor(hblup$u[valid], phenotype[valid], use="complete.obs")
  }
}

# Graphical representation of the results
###############################

# Boxplots

par(mfrow=c(1, 2), mar=c(8, 4, 2, 2))
# Predictive abilities obtained with GBLUP (genomic predictions)
boxplot(AccuGBLUP, ylab="Predictive abilities",
        ylim=c(-0.2, 1), las=2, col=c("blue", rep("lightblue", 7)),
        main="GBLUP (genomic prediction)", cex.axis=1)
abline(v=1.5)

# Predictive abilities obtained with GBLUP (genomic predictions)
boxplot(AccuHBLUP, ylab="Predictive abilities",
        ylim=c(-0.2, 1), las=2, col=c("red", rep("indianred", 7)),
        main="HBLUP (phenomic prediction)", cex.axis=1)
abline(v=1.5)

# Scatter plot HBLUP vs GBLUP

par(mfrow=c(1, 1), mar=c(4, 4, 4, 4))
plot(colMeans(AccuGBLUP), colMeans(AccuHBLUP),
     xlim=c(0, 1), ylim=c(0, 1),
     xlab="Predictive ability GBLUP (Genomic selection)",
     ylab="Predictive ability HBLUP (Phenomic selection)")
abline(a=0,b=1)
points(colMeans(AccuGBLUP)[1], colMeans(AccuHBLUP)[1],
       pch=22, col="red", bg="red") # highlight the reference environment
legend("bottomright", col=c("red", "black"),
       legend=c("Reference environment (with NIRS)","Other environments (without NIRS)"),
       pch=c(22, 1), cex=1)


# end of script