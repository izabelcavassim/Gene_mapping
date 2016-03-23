source("http://www.bioconductor.org/biocLite.R")
biocLite("multtest")
install.packages("gplots")
install.packages("LDheatmap")
install.packages("genetics")
install.packages("EMMREML")
install.packages("scatterplot3d")

library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

#rm(list=ls())
library(nlme)
setwd('C:/Users/MariaIzabel/Desktop/MASTER/GENE MAPPING/PROJECT/Project_data/Chrm2_test')

## Read marker data
geno<-read.csv('chr2.txt',header=F,sep='')
nAnim <- nrow(geno)
nSNP <- (ncol(geno)-1)/2

## Coding the genotypes 11=0, 12/21=1, 22=2
snpdose <- cbind(geno[,1],(geno[,2*(1:nSNP)]+geno[,2*(1:nSNP)+1])-2)
colnames(snpdose)[1] <- "ID"
rm(geno)

## Read phenotype file
pheno<-read.csv('phenotype.txt',sep='')

## Read pedigree data -- not necessary with exam data which has already merged pheno and pedigree
#pedigree <- read.table(file="dat1.ped", header=F, sep="\t", col.names=c("ID","sire","dam"))
#pheno <- merge(pedigree, pheno, by="ID", all.y=TRUE)
#head(pheno)

## Remove Father and Mother IDs from phenotype table (to fit specs from GAPIT)
pedigree = pheno[,1:3]
pheno = pheno[,-c(2:3)]

## Read position file
position<-read.csv('chr2.map',header=F,sep='')

## Rename the columns in position (to be the GM file)
colnames(position) <- c("Chr", "SNP", "BP")

## Rearrange columns in position to "SNP - Chr - BP"

# reorder by column name
#data <- data[c("A", "B", "C")]

#reorder by column index
position<- position[c(2,1,3)]
ID = snpdose[,1]
snpdose = snpdose[,-1] #had to remove ID from table to allow renaming the markers

## Rename the columns in snpdose (to be the GD file)
colnames(snpdose) = position[,1]
snpdose = cbind(ID, snpdose)

myGAPIT = GAPIT(Y = pheno, GD = snpdose, GM= position, PCA.total=3, Geno.View.output=FALSE,  group.from=0, group.to=2000,group.by=10)
