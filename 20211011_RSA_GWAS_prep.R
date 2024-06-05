library("reshape2")
library("data.table")

# load phenotype information
Pheno <- read.csv("/home/magdaj/Magdas_Pimp_SNPs/PimpsSNPS/final_SNPS/RSA_results/Pimp_pareto_RSA_for_GWAS.csv")
colnames(Pheno)

# load Genotype information
myG <- fread("/home/magdaj/Magdas_Pimp_SNPs/PimpsSNPS/final_SNPS/Pimp_01.HM.hmp.txt", sep="\t", header = TRUE)
head(myG)
dim(myG)


# issues with the genotypes not matching - GAPIT was previously quite fussy about it. So let's solve it
pheno_names <- unique(Pheno$genotype)
length(colnames(myG))
colnames(myG)
geno_names <- colnames(myG)[12:253]

pheno_names %in% geno_names
new_pheno <- subset(Pheno, Pheno$genotype %in% geno_names)
col.num <- colnames(myG) %in% pheno_names
col.num
NewGeno <- subset(myG, select = col.num)
head(NewGeno)
# let's double check if all is TRUE now?
colnames(NewGeno) %in% new_pheno$Taxa
# seems to work now fine

# let's SAVE it with the other part of HapMap file
myG_ess <- myG[,c(1:11)]
head(myG_ess)
myG1 <- cbind(myG_ess, NewGeno)
head(myG1)
dim(myG1)

# Let's add another chromosome:
# load Genotype information
myG <- fread("/home/magdaj/Magdas_Pimp_SNPs/PimpsSNPS/final_SNPS/Pimp_02.HM.hmp.txt", sep="\t", header = TRUE)
head(myG)
dim(myG)


# issues with the genotypes not matching - GAPIT was previously quite fussy about it. So let's solve it
geno_names <- colnames(myG)[12:253]

pheno_names %in% geno_names
col.num <- colnames(myG) %in% pheno_names
col.num
NewGeno <- subset(myG, select = col.num)
head(NewGeno)
# let's double check if all is TRUE now?
colnames(NewGeno) %in% new_pheno$Taxa
# seems to work now fine

# let's SAVE it with the other part of HapMap file
myG_ess <- myG[,c(1:11)]
head(myG_ess)
myG2 <- cbind(myG_ess, NewGeno)
head(myG2)
dim(myG2)

# Fuse Chromosome 1 and Chromosome 2 into one data:
myGu <- rbind(myG1, myG2)
