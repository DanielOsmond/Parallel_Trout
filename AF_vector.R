install.packages("devtools")
library(devtools)
remotes::install_github("petrelharp/templater")
install.packages("data.table","tidyverse","vcfR")
install.packages("magrittr")
install.packages("readr")
install.packages("dplyr")
install.packages("vctrs")
library("readr")
library("magrittr")
library("dplyr")
library("tidyverse")
library("vcfR")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("qvalue")

install.packages("remotes",repos = "http://cran.us.r-project.org")
remotes::install_github("JimWhiting91/afvaper")
library(afvaper,verbose = F)


lib = c("parallel","data.table","tidyr","dplyr","ggplot2","lostruct","ggpubr","vcfR", "gridExtra", "extrafont","devtools")
lapply(lib, library, character.only=TRUE)

setwd("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/mafs")

### read in global maf data per locus
chr_high_conf <- read_tsv("chr_mafs.tsv")
summary(chr_high_conf)

### can see that some loci have maf <0.05 and >0.95, we need to remove those

chr_high_conf %<>% filter(knownEM >= 0.05)
chr_high_conf %<>% filter(knownEM <= 0.95)
summary(chr_high_conf)

### it works (look at knownEM row) - retain ~4.9m sites / 5.6m 
pos1 <- subset(chr_high_conf, select = c(1,2))
write_tsv(pos1, "pos1.tsv")

### read in the per pop .maf files
### rename the column 'knownEM' with the pop name in nano/notepad before this stage 

AGA <- read_tsv("AGA.mafs")
AGB <- read_tsv("AGB.mafs")
AVB <- read_tsv("AVB.mafs")
AVD <- read_tsv("AVD.mafs")
CAM <- read_tsv("CAM.mafs")
CWA <- read_tsv("CWA.mafs")
CWB <- read_tsv("CWB.mafs")
CWC <- read_tsv("CWC.mafs")
HAY1 <- read_tsv("HAY1.mafs")
NEV2 <- read_tsv("NEV2.mafs")
OWE <- read_tsv("OWE.mafs")
TYA <- read_tsv("TYA.mafs")
TYD <- read_tsv("TYD.mafs")
TYE <- read_tsv("TYE.mafs")
WEA <- read_tsv("WEA.mafs")
WEF <- read_tsv("WEF.mafs")
YHA <- read_tsv("YHA.mafs")
YHE <- read_tsv("YHE.mafs")


#read in the good maf >0.05 snp loci 
pos1 <-read_tsv("pos1.tsv")

#create fields for position to match between dataframes
pos1$ID = paste0(pos1$chr, "__", pos1$pos)

AGA$ID = paste0(AGA$chromo, "__", AGA$position)
AGB$ID = paste0(AGB$chromo, "__", AGB$position)
AVB$ID = paste0(AVB$chromo, "__", AVB$position)
AVD$ID = paste0(AVD$chromo, "__", AVD$position)
CAM$ID = paste0(CAM$chromo, "__", CAM$position)
CWA$ID = paste0(CWA$chromo, "__", CWA$position)
CWB$ID = paste0(CWB$chromo, "__", CWB$position)
CWC$ID = paste0(CWC$chromo, "__", CWC$position)
#CWD$ID = paste0(CWD$chromo, "__", CWD$position)
HAY1$ID = paste0(HAY1$chromo, "__", HAY1$position)
#HAY2$ID = paste0(HAY2$chromo, "__", HAY2$position)
NEV2$ID = paste0(NEV2$chromo, "__", NEV2$position)
OWE$ID = paste0(OWE$chromo, "__", OWE$position)
TYA$ID = paste0(TYA$chromo, "__", TYA$position)
TYD$ID = paste0(TYD$chromo, "__", TYD$position)
TYE$ID = paste0(TYE$chromo, "__", TYE$position)
WEA$ID = paste0(WEA$chromo, "__", WEA$position)
WEF$ID = paste0(WEF$chromo, "__", WEF$position)
YHA$ID = paste0(YHA$chromo, "__", YHA$position)
YHE$ID = paste0(YHE$chromo, "__", YHE$position)

### tidy up columns unwanted - check in the dataframe to retain the pop maf and $ID

AGA_1 <- subset(AGA, select = -c(1,2,3,4))
AGB_1 <- subset(AGB, select = -c(1,2,3,4))
AVB_1 <- subset(AVB, select = -c(1,2,3,4))
AVD_1 <- subset(AVD, select = -c(1,2,3,4))
CAM_1 <- subset(CAM, select = -c(1,2,3,4))
CWA_1 <- subset(CWA, select = -c(1,2,3,4))
CWB_1 <- subset(CWB, select = -c(1,2,3,4))
CWC_1 <- subset(CWC, select = -c(1,2))
HAY1_1 <- subset(HAY1, select =-c(1,2))
NEV2_1 <- subset(NEV2, select =-c(1,2,3,4))
OWE_1 <- subset(OWE, select =-c(1,2,3,4))
TYA_1 <- subset(TYA, select =-c(1,2,3,4))
TYD_1 <- subset(TYD, select =-c(1,2,3,4))
TYE_1 <- subset(TYE, select =-c(1,2,3,4))
WEA_1 <- subset(WEA, select =-c(1,2,3,4))
WEF_1 <- subset(WEF, select =-c(1,2,3,4))
YHA_1 <- subset(YHA, select =-c(1,2,3,4))
YHE_1 <- subset(YHE, select =-c(1,2,3,4))

## merge loci locations and maf per pop ## 

new.table <- merge(pos1, AGA_1, by = "ID")
new.table %<>% merge(AGB_1, by = "ID")
new.table %<>% merge(AVB_1, by = "ID")
new.table %<>% merge(AVD_1, by = "ID")
new.table %<>% merge(CAM_1, by = "ID")
new.table %<>% merge(CWA_1, by = "ID")
new.table %<>% merge(CWB_1, by = "ID")
new.table %<>% merge(CWC_1, by = "ID")
new.table %<>% merge(HAY1_1, by = "ID")
new.table %<>% merge(NEV2_1, by = "ID")
new.table %<>% merge(OWE_1, by = "ID")
new.table %<>% merge(TYA_1, by = "ID")
new.table %<>% merge(TYD_1, by = "ID")
new.table %<>% merge(TYE_1, by = "ID")
new.table %<>% merge(WEA_1, by = "ID")
new.table %<>% merge(WEF_1, by = "ID")
new.table %<>% merge(YHA_1, by = "ID")
new.table %<>% merge(YHE_1, by = "ID")

AF_vector <- subset(new.table, select = -c(1))
AF_vector %<>% arrange(chr,pos)
write_tsv(AF_vector, "./AF_vector.tsv")

AF<-read_tsv("AF_vector.tsv")

### subsets ###
#1 = chr 
#2 = pos
#3 = AGA 
#4 = AGB 
#5 = AVB
#6 = AVD
#7 = CAM
#8 = CWA
#9 = CWB
#10 =NEV2
#11 = OWE
#12 = TYA
#13 = TYD 
#14 = TYE
#15 = WEA
#16 = WEF
#17 = YHA
#18 = YHE
#19=CWC
#20=HAY1
### change subset values based upon the included pops in the vector list ###


### AF-vapeR Running ###

### vector lists ###

high_select_list <- list(Anglesey = c("AGA","AGB"), # AGA = Anglesey metal, AGB = Anglesey clean
                         
                     #    Irish = c("AVB","OWE"), # OWE = Irish no metal, sep catchment. (I picked this one because no metal is probably better than less metal)
                         
                         Hayle = c("HAY1","CAM"), # HAY1 = Filthy Hayle (cornwall), CAM = Cornish clean. (Removed one of the duplicates, but either one is fine)
                         
                         Mid_Wales = c("YHE", "NEV2"), # NEV = Mid_wales clean
                         
                      #   N_Wales = c("CWA", "CWC"), # CWA = N-Wales metal, CWC = N-Wales lower metal (Pick whichever CWA one you think is best regional metal-clean pair)
                         
                         NE_England = c("TYE", "TYA"), # TYA = NE-England clean, TYD = NE-England metal
                         
                         NE_England2 = c("WEA", "WEF")) #WEA = NE-England metal, WEF = NE clean. (Consider maybe removing this pair if it is very close to other NE ENgland one?)

vector_list_UK <- list(Anglesey = c("AGA","AGB"), # AGA = Anglesey metal, AGB = Anglesey clean
                     
                     Irish = c("AVB","OWE"), # OWE = Irish no metal, sep catchment. (I picked this one because no metal is probably better than less metal)
                     
                     Hayle = c("HAY1","CAM"), # HAY1 = Filthy Hayle (cornwall), CAM = Cornish clean. (Removed one of the duplicates, but either one is fine)
                     
                     Mid_Wales = c("YHE", "NEV2"), # NEV = Mid_wales clean
                     
                     N_Wales = c("CWA", "CWC"), # CWA = N-Wales metal, CWC = N-Wales lower metal (Pick whichever CWA one you think is best regional metal-clean pair)
                     
                     NE_England = c("TYD", "TYA"), # TYA = NE-England clean, TYD = NE-England metal
                     
                     NE_England2 = c("WEA", "WEF")) #WEA = NE-England metal, WEF = NE clean. (Consider maybe removing this pair if it is very close to other NE ENgland one?)

vector_list_NE <- list(TY1 = c("TYD", "TYA"),
                       TY2 = c("TYE", "TYA"),
                       WE1 = c("WEA", "WEF"))

vector_list_cymru <- list(Anglesey = c("AGA", "AGB"),
                          Mid_Wales = c("YHE", "NEV2"), #which YH site to use??
                          N_wales = c("CWA", "CWC"),
                          N_wales_2 =c("CWB", "CWC"))
                          

vector_list_Cu <- list(Anglesey = c("AGA", "AGB"),
                       Irish = c("AVB","OWE"),
                       Hayle = c("HAY1","CAM"))

vector_list_Pb<- list(Mid_Wales = c("YHE", "NEV2"),
                      TY2 = c("TYE", "TYA"),
                      N_wales_2 =c("CWB", "CWC"),
                      WE1 = c("WEA", "WEF"))
  
  ### calculate AF_vectors ### 

window_snps=200
chrs<-as.list(unique(AF$chr))

AF_vector_57 <- data.frame(subset(AF, chr =="NC_042957.1"))
AF_vector_58 <- data.frame(subset(AF, chr =="NC_042958.1"))
AF_vector_59 <- data.frame(subset(AF, chr =="NC_042959.1"))
AF_vector_60 <- data.frame(subset(AF, chr =="NC_042960.1"))
AF_vector_61 <- data.frame(subset(AF, chr =="NC_042961.1"))
AF_vector_62 <- data.frame(subset(AF, chr =="NC_042962.1"))
AF_vector_63 <- data.frame(subset(AF, chr =="NC_042963.1"))
AF_vector_64 <- data.frame(subset(AF, chr =="NC_042964.1"))
AF_vector_65 <- data.frame(subset(AF, chr =="NC_042965.1"))
AF_vector_66 <- data.frame(subset(AF, chr =="NC_042966.1"))
AF_vector_67 <- data.frame(subset(AF, chr =="NC_042967.1"))
AF_vector_68 <- data.frame(subset(AF, chr =="NC_042968.1"))
AF_vector_69 <- data.frame(subset(AF, chr =="NC_042969.1"))
AF_vector_70 <- data.frame(subset(AF, chr =="NC_042970.1"))
AF_vector_71 <- data.frame(subset(AF, chr =="NC_042971.1"))
AF_vector_72 <- data.frame(subset(AF, chr =="NC_042972.1"))
AF_vector_73 <- data.frame(subset(AF, chr =="NC_042973.1"))
AF_vector_74 <- data.frame(subset(AF, chr =="NC_042974.1"))
AF_vector_75 <- data.frame(subset(AF, chr =="NC_042975.1"))
AF_vector_76 <- data.frame(subset(AF, chr =="NC_042976.1"))
AF_vector_77 <- data.frame(subset(AF, chr =="NC_042977.1"))
AF_vector_78 <- data.frame(subset(AF, chr =="NC_042978.1"))
AF_vector_79 <- data.frame(subset(AF, chr =="NC_042979.1"))
AF_vector_80 <- data.frame(subset(AF, chr =="NC_042980.1"))
AF_vector_81 <- data.frame(subset(AF, chr =="NC_042981.1"))
AF_vector_82 <- data.frame(subset(AF, chr =="NC_042982.1"))
AF_vector_83 <- data.frame(subset(AF, chr =="NC_042983.1"))
AF_vector_84 <- data.frame(subset(AF, chr =="NC_042984.1"))
AF_vector_85 <- data.frame(subset(AF, chr =="NC_042985.1"))
AF_vector_86 <- data.frame(subset(AF, chr =="NC_042986.1"))
AF_vector_87 <- data.frame(subset(AF, chr =="NC_042987.1"))
AF_vector_88 <- data.frame(subset(AF, chr =="NC_042988.1"))
AF_vector_89 <- data.frame(subset(AF, chr =="NC_042989.1"))
AF_vector_90 <- data.frame(subset(AF, chr =="NC_042990.1"))
AF_vector_91 <- data.frame(subset(AF, chr =="NC_042991.1"))
AF_vector_92 <- data.frame(subset(AF, chr =="NC_042992.1"))
AF_vector_93 <- data.frame(subset(AF, chr =="NC_042993.1"))
AF_vector_94 <- data.frame(subset(AF, chr =="NC_042994.1"))
AF_vector_95 <- data.frame(subset(AF, chr =="NC_042995.1"))
AF_vector_96 <- data.frame(subset(AF, chr =="NC_042996.1"))


vector_list <- list(AF_vector_57,AF_vector_58,AF_vector_59,AF_vector_60,
                    AF_vector_61,AF_vector_62,AF_vector_63,AF_vector_64,
                    AF_vector_65,AF_vector_66,AF_vector_67,AF_vector_68,
                    AF_vector_69,AF_vector_70,AF_vector_71,AF_vector_72,
                    AF_vector_73,AF_vector_74,AF_vector_75,AF_vector_76,
                    AF_vector_77,AF_vector_78,AF_vector_79,AF_vector_80,
                    AF_vector_81,AF_vector_82,AF_vector_83,AF_vector_84,
                    AF_vector_85,AF_vector_86,AF_vector_87,AF_vector_88,
                    AF_vector_89,AF_vector_90,AF_vector_91,AF_vector_92,
                    AF_vector_93,AF_vector_94,AF_vector_95,AF_vector_96)


### loop through the vector list of individual chromosomes and calc AF vectors ###


all_chr_res <- lapply(1:40,function(i){
  
  ### AF_subsets (select 1&2 for chr, pos, then pops based on index above)
  
  #UK_AF <- subset(vector_list[[i]], select = -c(14,9,6)) #UK
  #AF <- subset(vector_list[[i]], select = c(1,2,12,13,14,15,16)) #NE
  #AF <- subset(vector_list[[i]], select = c(1,2,3,4,8,9,10,19,18)) #cymru
  #AF <- subset(vector_list[[i]], select = c(1,2,3,4,5,11,20,7)) #Cu
  #AF <- subset(vector_list[[i]], select = c(1,2,3,4,20,7,18,10,14,12,15,16))
  AF <- subset(vector_list[[i]], select = c(1,2,18,10,14,12,9,19,15,16))
  ### issue around $pos being character string rather than numeric - fix below ### 
  
  AF$pos=as.integer(AF$pos)
  AF=AF[order(AF$pos),]
  
  chr_AF_input <- calc_AF_vectors(vcf = AF,
                                  window_size = window_snps,
                                  vectors = vector_list_Pb,
                                  n_cores = 1,
                                  data_type = "freq")
  
  # How many permutations to run
  null_perm_N = 1000
  
  # Calculate Allele Frequency Change Vector Matrices
  chr_null_input <- calc_AF_vectors(vcf = AF,
                                window_size = window_snps,
                                vectors = vector_list_Pb,
                                n_cores = 1,
                                null_perms = null_perm_N,
                                data_type = "freq")
  
  ## We could save these to some temporary file, e.g.
  saveRDS(list(chr_AF_input,chr_null_input),paste0("Pb","chr",i,"_AFV.rds"))
  
  # Return our results
  return(list(chr_AF_input,chr_null_input))
  
  
})

## If saved as .rds, we can read these back in and merge to single lists
all_chr_res <- lapply(1:40,function(i) readRDS(paste0("NE","chr",i,"_AFV.rds")))

# To fetch all of  chr AFV, we take the first element of each list element
# Note: the merge_eigen_res() func is the same as unlist(,recursive=F)
AF_input <- merge_eigen_res(lapply(all_chr_res,'[[',1))

# All null, we take the second element of each list element
null_input <- merge_eigen_res(lapply(all_chr_res,'[[',2))

# We now have our whole (3 chr) genome's worth of AFV matrices and null matrices in a single input
c(head(names(AF_input)),tail(names(AF_input)))

# Perform eigen analysis
eigen_res <- lapply(AF_input,eigen_analyse_vectors)

# View chromosomal regions:
head(names(eigen_res))

# View eigenvalue distribution of first matrix
eigen_res[[1]]$eigenvals

# View eigenvector loadings of first matrix
eigen_res[[1]]$eigenvecs

# View head of SNP scores
head(eigen_res[[1]]$A_matrix)

# Get cutoffs for 95%, 99% and 99.9%
null_cutoffs <- find_null_cutoff(null_input,cutoffs = c(0.95,0.99,0.999))
null_cutoffs

# Calculate p-vals
pvals <- eigen_pvals(eigen_res,null_input)

# Showpvals
head(pvals)

# Plot the raw eigenvalues, and visualise the cutoff of 99%
all_plots <- eigenval_plot(eigen_res,cutoffs = null_cutoffs[,"99%"])

# Show the plots for eigenvalue 1
all_plots[[1]] + theme(axis.text.x=element_text(size=2))
ggsave("NE_eig1.png", width = 30, height =4, dpi =300, bg ="white")

# Plot empirical p-values, -log10(p) of 2 ~ p=0.01, 3 ~ p=0.001 etc.
all_plots_p <- eigenval_plot(eigen_res,null_vectors = null_input,plot.pvalues = T)
all_plots_p[[3]]

library(ggplot2,verbose = F)

# Pull the figure 
eig1_pval_fig <- all_plots_p[[1]]

# Edit
eig1_pval_fig + theme(title = element_blank()) + geom_step(colour="red2")

# Recall the use of find_null_cutoffs() to fetch a matrix of cutoffs...
# null_cutoffs

# Find significant windows above 99.9% null permutation
significant_windows <- signif_eigen_windows(eigen_res,null_cutoffs[,"99.9%"])

# Display 'outliers'
significant_windows

# Summarise parallel evolution in windows that are significant on eigenvector 1
eig1_parallel <- summarise_window_parallelism(window_id = significant_windows[[1]],
                                              eigen_res = eigen_res,
                                              loading_cutoff = 0.3,
                                              eigenvector = 1)


# Show results
head(eig1_parallel)

options(max.print = 1000)

### save the output ###
### stretch the console to get full width ### 
### RENAME ### 
### RENAME ### 
### RENAME ###
sink("NE_eig1_99.9%.txt")
print(eig1_parallel)
sink()

write.table(eig1_parallel, file="Cu_eig1_99.9%.txt", row.names = T, col.names = T)

### second eigenvector ###
# Summarise
eig2_parallel <- summarise_window_parallelism(window_id = significant_windows[[2]],
                                              eigen_res = eigen_res,
                                              loading_cutoff = 0.3,
                                              eigenvector = 2)

# Show results
head(eig2_parallel)



### save the output ### 
### RENAME ### 
### RENAME ### 
### RENAME ### 
sink("Cu_eig2_99%.txt")
print(eig2_parallel)
sink()

write.table(eig2_parallel, file="cymru_eig2_99%.txt", row.names = T , col.names = T)

## explore candidate regions ## 
## This a matrix with a row per SNP and a column per eigenvector ##
A_mat <- eigen_res[[1]]$A_matrix
head(A_mat)

### These can be plotted by separating out the chr and pos values stored in rownames(A_mat)###
### Note that here we have not polarised the allele frequencies in any way, 
### by default the analysis just tracks the REF allele as labelled in the VCF, 
### so both positive and negative scores of association are of interest. 
### Here then we can just plot the absolute score:

to_plot <- data.frame(snp=rownames(A_mat),
                      eig1_score=A_mat[,1])

to_plot <- to_plot %>% tidyr::separate("snp",into=c("NA", "chr","pos"),sep="_")
to_plot$pos <- as.integer(to_plot$pos)

eig1_pos<-ggplot(to_plot,aes(x=pos,y=abs(eig1_score)))+
  geom_point()+
  labs(y="Eig1 Score",x="Pos (bp)")

eig1_pos# + xlim(0,12000000)
ggsave("UK_eig1_pos_plot.png", dpi=300, width=10, height=6)

