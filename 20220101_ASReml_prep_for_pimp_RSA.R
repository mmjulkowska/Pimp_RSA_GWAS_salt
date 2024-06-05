setwd("/home/magdaj/Magdas_Pimp_SNPs/PimpsSNPS/final_SNPS/")
list.files(pattern=".csv")
Geno <- read.csv("Pimp_all_chromosomes.csv", header = TRUE)
head(Geno)
colnames(Geno) <- Geno[1,]
dim(Geno)
Geno2 <- Geno[2:708546,2:205]
head(Geno2)
Pheno <- read.csv("RSA_results/Pimp_pareto_RSA_for_GWAS.csv")
Pheno

pheno_names <- unique(Pheno$genotype)
length(pheno_names)
length(colnames(Geno2))
colnames(Geno2)
geno_names <- colnames(Geno2)[12:204]

pheno_names %in% geno_names
geno_names %in% pheno_names
col.num <- colnames(Geno2) %in% pheno_names
col.num
NewGeno <- subset(Geno2, select = col.num)
head(NewGeno)
dim(NewGeno)


# let's double check if all is TRUE now?
colnames(NewGeno) %in% Pheno$genotype
# seems to work now fine

# Let's save it with another part of the HapMap file:
myG_ess <- Geno2[,c(1:11)]
head(myG_ess)
myGeno <- cbind(myG_ess, NewGeno)
head(myGeno)
dim(myGeno)

# looks good to me! 
# Let's try a quick GAPIT if it all works well:
colnames(Pheno)
dim(Pheno)
head(Pheno)
Pheno2 <- Pheno[,2:538]
head(Pheno2)
head(myGeno)
unique(myGeno$chrom)

library(GAPIT3)
myGAPIT <- GAPIT(
  Y=Pheno[,c(1,2)], #fist column is individual ID, the third columns is days to pollination
  G=myGeno,
  PCA.total=3)


# OK - since we need to get the data in a right format - 
# let's create a file where we will put in the numeric values
dim(myGeno)
head(myGeno)
myGugu <- myGeno[,12:193]
MyGeno2 <- myGugu
colnames(MyGeno2)

temp.list <- as.numeric(as.factor(myGugu[1,]))
MyGeno2[1,] <- gsub(1, 0, temp.list)
head(MyGeno2)

# We also need a file with mapAdaptom.txt where we have SNP name + Chr + Pos information
mapAdaptom <- myGeno[,c(1,3,4)]
colnames(mapAdaptom)[1] <- "SNP"
head(mapAdaptom)
tail(mapAdaptom)
write.table(mapAdaptom, "mapAdaptom.txt", row.names = FALSE)
unique(mapAdaptom$chrom)
dim(MyGeno2)

# Now - let's loop it!
for(i in 1:708545){
  temp.list <- as.numeric(as.factor(myGugu[i,]))
  MyGeno2[i,] <- gsub(1, 0, temp.list)
}

head(MyGeno2)


write.csv(MyGeno2, "RSA_GenotypeNumeric.csv", row.names = F)

#  Now I should go on and execute the following:

head(Pheno2)

colnames(Pheno)[1] <- "ecot_id"
Pheno$ecot_id %in% colnames(MyGeno2)
my_YY <- Pheno
my_YY$ecot_id <- gsub("M", "", my_YY$ecot_id)
my_YY$ecot_id %in% colnames(MyGeno2)

my_X <- MyGeno2
colnames(my_X) <- gsub("M", "", colnames(my_X))
my_YY

my_YY$ecot_id %in% colnames(my_X)

# I am still missing my kinship - but let's use the old kinship file:
K <- read.csv("GAPIT.Kin.VanRaden.csv",header=F)
dim(K)
K1 <- K[,c(2:183)]
K1

length(colnames(my_X))
dim(K1)
colnames(K1) <- colnames(my_X)
rownames(K1) <- colnames(my_X)
my_K <- as.matrix(K1)

my_SNP <- read.table("mapAdaptom.txt", header = T)
head(my_SNP)
tail(my_SNP)
dim(my_SNP)
dim(my_X)

rownames(my_X) <- my_SNP$SNP

my_Xt <- t(my_X)
head(my_Xt)
dim(my_Xt)
my_Xt[1:5,1:5]
# ok - I have to change it to numeric
X_num <- matrix(as.numeric(my_Xt),    # Convert to numeric matrix
                ncol = ncol(my_Xt))
X_num[1:5,1:5]
colnames(X_num) <- colnames(my_Xt)
rownames(X_num) <- as.numeric(rownames(my_Xt))
my_YY$ecot_id <- as.numeric(my_YY$ecot_id)


load("Pimp.TPA.GWAS.data.RData")
library(asreml)
source("t250_gwas.R")
head(Y)

head(my_YY)
dim(my_YY)

dim(my_X)
dim(my_SNP)
unique(my_SNP$chrom)
Y <- my_YY
K <- my_K
X <- X_num

my_YY$ecot_id %in% colnames(K)
my_YY$ecot_id %in% rownames(X)
dim(my_YY)
dim(X)

new_YY <- subset(my_YY, my_YY$ecot_id %in% rownames(X))
dim(new_YY)
head(K)
new_YY$ecot_id %in% colnames(K)
new_YY$ecot_id %in% rownames(X)

new_YY[1:5, 1:5]
X[1:5, 1:5]

X <- matrix(as.numeric(X),    # Convert to numeric matrix
                ncol = ncol(X))
X[1:5,1:5]
colnames(X) <- colnames(my_Xt)
rownames(X) <- as.numeric(rownames(my_Xt))
dim(new_YY)

new_YY2 <- new_YY[,c(1:12,33:42,63:78,99:108,129:144,165:174,195:210,231:240,261:276,297:306,
                     327:342,363:372,393:408,429:438,459:474, 495:504, 525:538)]

colnames(new_YY2)
# Let's try if this is going to work
for(i in 2:218){
  t2029_gwas(new_YY2,n=i)  
}  


# After we have done this - we still need to generate file that has the SNP frequencies and such. 
# For this purpose we need to run this:
names <- c(text="SNP", "MAC", "MAF")
SNPskies <- data.frame()
SNPskies
for (k in names) SNPskies[[k]] <- as.character()

SNPskies[1,1] <- colnames(X)[1]
SNPskies[1,2] <- sum(X[,1])/2
SNPskies[1,3] <- as.numeric(SNPskies[1,2])/220
head(SNPskies)

# let's loop it:
dim(X)
# ADJUST HERE SNP numbers
for(i in 2:708545){
  SNPskies[i,1] <- colnames(X)[i]
  MAC_temp <- sum(X[,i])/2
  if(as.numeric(MAC_temp) > 110){
    SNPskies[i,2] <- 222 - MAC_temp
    SNPskies[i,3] <- 1 - as.numeric(SNPskies[i,2])/222
  } else {
    SNPskies[i,2] <- as.numeric(MAC_temp)
    SNPskies[i,3] <- as.numeric(SNPskies[i,2])/222
  }
}

write.csv(SNPskies, "SNP_freq_RSA.csv", row.names = FALSE)

save(new_YY2, my_K, my_SNP, X_num, file = "Pimp.RSA.GWAS.data.RData")


# OK - now let's try to re-run this with calculations for effect size
# in order for this to work - I had to restart R and load all the data again from the beginning
load("Pimp.RSA.GWAS.data.RData")
library(asreml)
source("t250_gwas.R")
head(new_YY2)
colnames(new_YY2)
Y <- new_YY2[,c(1:6, 10:13, 17:31, 36:39, 43:57, 62:65, 69:83, 88:91, 95:109, 114:117, 121:135, 140:161, 166:187, 192:195, 199:218)]
X_num[1:5,1:5]
X <- X_num
K <- my_K
dim(Y)

ncol(Y)
i=2
t2029_gwas(Y,n=i)  

for(i in 3:169){
  t2029_gwas(Y,n=i)  
}