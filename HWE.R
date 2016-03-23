setwd('C:/Users/MariaIzabel/Desktop/MASTER/GENE MAPPING/PROJECT/Project_data')

each_chromosome = vector()
names_not_HWE = vector()
chromosomes = c('chr1.txt', 'chr2.txt', 'chr3.txt', 'chr4.txt', 'chr5.txt')
physical_map = c('chr1.map', 'chr2.map', 'chr3.map', 'chr4.map', 'chr5.map')

for(k in chromosomes){
  for(x in physical_map){
  geno = read.table(k)
  snp_names = read.table(x, as.is=TRUE)
  snp_names <- snp_names[,2]

  nSNP <- (ncol(geno)-1)/2
  ## Coding the genotypes 11=0, 12/21=1, 22=2
  snpdose <- cbind(geno[,1],(geno[,2*(1:nSNP)]+geno[,2*(1:nSNP)+1])-2)
  colnames(snpdose)[1] <- "ID"
  rm(geno)

  results = vector()
  names = vector()
  colnames(snpdose) = seq(1,1998)
  for(i in 2:ncol(snpdose)){
    N0<-length(snpdose[,i][snpdose[,i] == 0])
    N1<-length(snpdose[,i][snpdose[,i] == 1])
    N2<-length(snpdose[,i][snpdose[,i] == 2])
    total = N0 + N1 + N2
    p = mean(snpdose[,2])/2
    q = 1 - p
    HWe_expected_AA = p^2*total
    HWe_expected_Aa = 2*p*q*total
    HWe_expected_aa = q^2*total
    observed = c(N0, N1, N2)
    expected = c(HWe_expected_AA, HWe_expected_Aa, HWe_expected_aa)
  
    test = cbind(observed, expected)
    Chi2Sum = sum(((test[,1] - test[,2])^2)/test[,2])
    Pvalue = 1-pchisq(q=Chi2Sum,df = 2)
    results[i] = Pvalue
  }

# The number of SNPs that are in HWE is:
each_chromosome[k] = sum(results <= 0.05, na.rm='TRUE')
if(k == 'chr1.txt'){
names_not_HWE = append(snp_names[which(results >= 0.05)], names_not_HWE)
}
  }
}

sort(names_not_HWE)
