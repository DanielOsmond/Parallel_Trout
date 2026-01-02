### 
# SFS Utilities 

setwd("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/sfs")
library(readr)

AGA_sfs<-scan(paste0("./AGA.sfs"))
barplot(AGA_sfs)
barplot(AGA_sfs[c(-1, -length(AGA_sfs))], ylim=c(0,600000))

CWC_sfs<-scan(paste0("./CWC.sfs"))
barplot(CWC_sfs)
barplot(AGB_sfs[c(-1, -length(AGB_sfs))], ylim =c(0,600000))


### evomics script ### 
##run in R
#plot the results
nnorm <- function(x) x/sum(x)
#expected number of sites with 1:20 derived alleles
#res <- rbind(
#  AGA=scan("AGA.sfs")[-1],
#  AGB=scan("AGB.sfs")[-1],
#  AVB=scan("AVB.sfs")[-1],
#  AVD=scan("AVD.sfs")[-1],
#  CAM=scan("CAM.sfs")[-1],
#  CWA=scan("CWA.sfs")[-1],
#  CWB=scan("CWB.sfs")[-1],
#  CWC=scan("CWC.sfs")[-1],
#  HAY1=scan("HAY1.sfs")[-1],
#  NEV2=scan("NEV2.sfs")[-1],
#  OWE=scan("OWE.sfs")[-1],
#  TYA=scan("TYA.sfs")[-1],
#  TYD=scan("TYD.sfs")[-1],
#  TYE=scan("TYE.sfs")[-1],
#  WEA=scan("WEA.sfs")[-1],
#  WEF=scan("WEF.sfs")[-1],
#  YHA=scan("YHA.sfs")[-1],
#  YHE=scan("YHE.sfs")[-1],
#)
colnames(res) <- 1:20

grouped <- read.table("grouped_sfs.csv", sep= ",", header = F)
grouped_mat <- as.matrix(grouped[,2:22])
names(grouped) <-c("pop", 0:20)

# density instead of expected counts
#res <- t(apply(res,1,nnorm))
grouped_freq <- subset(grouped, pop =="AGA" | pop== "AVB" | pop=="CAM" |pop== "HAY1" | pop=="TYA")
grouped_freq <- t(apply(grouped[,2:22],1,nnorm))

#plot the none ancestral sites
barplot(grouped_freq,beside=T,legend=c("AGA", "AVB", "CAM", "HAY1", "TYA"),names=1:20,main="realSFS non ancestral sites",args.legend = list(x="topright"))

#plot the polymorphic sites.
resPoly <- t(apply(grouped_freq[,-20],1,nnorm))
barplot(resPoly,beside=T,legend=c("AGA", "AVB", "CAM", "HAY1", "TYA"),names=1:19,main="realSFS polymorphic sites")

#due the very limited amount of sites
#downsample to 5 individuals (10 chromosome) and exclude fixed derived
downsampleSFS <- function(x,chr){ #x 1:2n , chr < 2n
  n<-length(x)
  mat <- sapply(1:chr,function(i) choose(1:n,i)*choose(n- (1:n),chr-i)/choose(n,chr))
  nnorm( as.vector(t(mat) %*% x)[-chr] )
}
resDown <- t(apply(grouped_mat,1,downsampleSFS,chr=10))
barplot(resDown,beside=T,legend=c("AGA",  "AGB",  "AVB",  "AVD",
"CAM",  "CWA",  "CWB",  "CWC", "HAY1", "NEV",  "OWE",  "TYA"),
names=1:9,main="realSFS downsampled polymorphic sites")

### calculate statistics ### 

AGA <- scan("AGA.sfs")
AGB <- scan("AGB.sfs")
AVB <- scan("AVB.sfs")
AVD <- scan("AVD.sfs")
CAM <- scan("CAM.sfs")
CWA <- scan("CWA.sfs")
CWB <- scan("CWB.sfs")
CWC <- scan("CWC.sfs")
HAY1 <- scan("HAY1.sfs")
NEV2 <- scan("NEV2.sfs")
OWE <- scan("OWE.sfs")
TYA <- scan("TYA.sfs")
TYD <- scan("TYD.sfs")
TYE <- scan("TYE.sfs")
WEA <- scan("WEA.sfs")
WEF <- scan("WEF.sfs")
YHA <- scan("YHA.sfs")
YHE <- scan("YHE.sfs")

pop <- HAY1

nSites<-sum(pop)   #Number of sites where we have data
nSeg<-sum(pop[c(-1,-21)])    #Number of segregating sites
an <- function(n) sum(1/1:(n-1))
thetaW <- nSeg/an(20) # Wattersons Theta
pop_n<-thetaW / 8e-9 / nSites / 4 # effective population size

HAY1_results<-as.data.frame.list(list(nSites, nSeg, thetaW,pop_n))
print(HAY1_results[1,])
