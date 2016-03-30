setwd('C:/Users/MariaIzabel/Desktop/MASTER/GENE MAPPING/PROJECT/Project_data')
chr1 = read.csv('GAPIT..phenotype.GWAS.Results_chr1.txt', header = TRUE, sep=',')
chr2 = read.csv('GAPIT..phenotype.GWAS.Results_chr2.csv', header = TRUE, sep=',')
chr3 = read.csv('GAPIT..phenotype.GWAS.Results_chr3.csv', header = TRUE, sep=',')
chr4 = read.csv('GAPIT..phenotype.GWAS.Results_chr4.csv', header = TRUE, sep=',')
chr4$Chromosome[chr4$Chromosome==2] <- as.numeric(4) #a little mistake :P
chr5 = read.csv('GAPIT..phenotype.GWAS.Results_chr5.csv', header = TRUE, sep=',')

results = rbind(chr1, chr2, chr3, chr4, chr5)
require(qqman)
head(results)
results <-(subset(results, select=c(SNP, Chromosome, Position, P.value)))
head(results)
colnames(results) <- c("SNP", "CHR", "BP", "P")      ## these are the names used in 'manhattan'

tail(results)
#How many SNPs on each chromosome?
as.data.frame(table(results$CHR))
#boring version
pdf('Manhattan_boring.pdf')
manhattan(results,  suggestiveline = F, genomewideline = F)

#cool colors:
rainbow(5)
pdf('Manhattan_colorful.pdf')
manhattan(results, main = "Manhattan Plot", ylim = c(0, 8), cex = 0.95, chrlabs= c('1','2','3','4', '5'), 
    col = c( "#FF0000FF", "#CCFF00FF", "#00FF66FF", "#0066FFFF", "#CC00FFFF"), suggestiveline = F, genomewideline = F )
dev.off()

#highligthing snps of interest:
SNPs_interesting_chr1 = c('chr1:2850000', 'chr1:2900000', 'chr1:2950000',
 'chr1:2700000', 'chr1:2650000', 'chr1:3300000', 'chr1:3500000', 'chr1:3200000', 'chr1:1550000')

pdf('manhattan_chrm1.pdf')
manhattan(subset(results, CHR == 1), highlight= SNPs_interesting_chr1, main = "Interesting SNPs")
dev.off()
# qqplot of the entire data




pdf('qqplot_entire_data.pdf')
qq(results$P, main = "Q-Q plot of GWAS p-values", xlim = c(0, 7), ylim = c(0, 
    12), pch = 20, col = "blue4", cex = 1.5, las = 1)

