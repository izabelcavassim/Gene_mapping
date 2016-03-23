## Program for single marker analysis using linear mixed model without fitting relationship  among individuals
   rm(list=ls())
   require('nlme')
   setwd('C:/Users/MariaIzabel/Desktop/MASTER/GENE MAPPING/PROJECT/Project_data')
 ## Read marker data
   geno <- read.table(file='chr1.txt',header=FALSE )
   nAnim <- nrow(geno)
   nSNP <- (ncol(geno)-1)/2

## Coding the genotypes 11=0, 12/21=1, 22=2
   snpdose <- cbind(geno[,1],(geno[,2*(1:nSNP)]+geno[,2*(1:nSNP)+1])-2)
   colnames(snpdose)[1] <- "ID"
   rm(geno)
## Read phenotype file
   pheno <- read.table(file='phenotype.txt', header=T, sep="\t")
## Read pedigree data
   #pedigree <- read.table(file="dat1.ped", header=F, sep="\t", col.names=c("ID","sire","dam"))
   #pheno <- merge(pedigree, pheno, by="ID", all.y=TRUE)
   head(pheno)
  
  #####################
 
  ### sire model analysis
       outTable <- matrix(NA,nrow=nSNP,ncol=6)

  ### Run Father model for all SNPs  
     for ( i in 1:nSNP) { 
        mydata <- data.frame(anim.ID=as.factor(pheno[,1]), father.ID=as.factor(pheno[,2]),
        phe=pheno[,4])
        mydata$mysnp <- snpdose[,i+1]
     
       if (length(levels(as.factor(mydata$mysnp))) > 1) {
###Model 
### Running father as fixed effect 

      # lmF <- lm( phe ~ 1 + sire.ID + mysnp, data=mydata)    
      # lmp <- anova(lmF)["mysnp","Pr(>F)"]
      # DF <- lmF$df.residual 
      # Mlogp <- -log10(lmp)          # Mlogp is -log10(p)
      # outTable[i,] <- c(j,i, coef(summary(lmF))["mysnp",], Mlogp)
        

### Running sire as random effect 
   
       fm <-lme(fixed=phe ~ 1 + mysnp, random=~1|father.ID, data=mydata, control=list(returnObject=TRUE))  ## if not converging
       pout <- summary(fm)$tTable[2,]   
       mlogp <- -log10(pout[5])
       outTable[i,] <- c(i,  pout[c(1,2, 4,5)], mlogp)
        }  else outTable[i,] <- c(i, rep(NA,5))                 
                  
       }     # for SNP

## Write output
   colnames(outTable) <- c("SNPno","beta","SE","t_value","p_value","minuslog10p")
   write.table(outTable,file="gwas-sireRandom.txt",append=F,quote=FALSE, row.name=FALSE, sep='\t')

  
  
##Plot GWAS results
source("manhattanPlot.r")      ## Function for manhattan plot


gwas <- read.table(file="gwas-sireRandom.txt", header=T)
map <- read.csv('chr1.map', sep='', header=F)

comb <- cbind(gwas,map)

   results <-(subset(comb, select=c(SNPno, V1, V3, p_value)))
   colnames(results) <- c('SNP', 'CHR', 'BP', 'P') 
   head(results)
   
  #postscript(file=paste("manhattan_CH4",tr,".eps",sep=""), width = 600, height = 300)
  pdf(file="manhattan-sireRandom.pdf")
   manhattan(results,pch=20, cex=0.3, main="Manhattan Plot", genomewideline=F, suggestiveline=F)
   abline(h=-log10(0.05/(nrow(results))), col="red")
  dev.off()

  ## QQ plot
  pdf(file="qq-sireRandom.pdf")
   qq(results$P)
  dev.off()
  
  # taking off the NA entries in our data 
  results = results[!is.na(results$P),]
  ###Get lambda value
  chisq <- qchisq(1-results$P,1)
  lamda <- median(chisq)/qchisq(0.5,1)
  print(lamda)