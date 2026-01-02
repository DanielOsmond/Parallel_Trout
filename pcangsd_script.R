install.packages("RcppCNPy")
install.packages("tidyverse")
install.packages("ggfortify")#
library("ggfortify")
library("RcppCNPy")
library("tidyverse")
library("cowplot")
library("MASS")
library("memisc")

setwd("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/pcangsd_plotting")

cov <- as.matrix(read.table("pcangsd_linkage_pruned_PCA.cov"))
samples <- read.csv("sample_info.csv", header = T)
ind <- as.matrix(samples$ind)
pop <- as.matrix(samples$pop)

mme.pca <- eigen(cov)

eigenvectors = mme.pca$vectors #extract eigenvectors 
pca.vectors = as_tibble(cbind(pop, data.frame(eigenvectors))) #combine with our population assignments

pca = ggplot(data = pca.vectors, aes(x=X1, y=X2, colour = pop)) + geom_point() +
  theme_classic()
plot(pca)
#ggsave("depth_unpruned_pcangsd_1v3.pdf", plot =pca, width =4, height =6)

pca.eigenval.sum = sum(mme.pca$values) #sum of eigenvalues
varPC1 <- (mme.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (mme.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (mme.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (mme.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4

source("individual_pca_functions.R")

pca_12 = PCA(cov, ind, pop, 1, 2, show.point=T,show.label= F, show.line=T) 
DAPC_12 = DAPC(n= 16, 1, 2, show.point = T, show.label = T, show.line = F)
pca_13 = PCA(cov, ind, pop, 1, 3, show.point=T,show.label= T, show.line=F)
DAPC_13 = DAPC(n= 15, 1, 3, show.point = T, show.label = T, show.line = F)

ggsave("linkage_pruned_pca_nolabs_1v2.pdf", plot =pca_12, width =6, height =6)
ggsave("linkage_pruned_pca_1v3.pdf", plot =pca_13, width =6, height =6)
ggsave("linkage_pruned_DAPC_1v2.pdf", plot =DAPC_12, width =6, height =6)
ggsave("linkage_pruned_DAPC_1v3.pdf", plot =DAPC_13, width =6, height =6)






