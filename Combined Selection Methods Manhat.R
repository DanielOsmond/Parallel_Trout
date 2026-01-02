##################################################
#### selection scan figure ###########
##################################################

setwd("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics")

library("readr")
library("magrittr")
library("dplyr")
library("tidyverse")
library("vcfR")
library(afvaper,verbose = F)
lib = c("parallel","data.table","tidyr","dplyr","ggplot2","ggpubr","vcfR", "gridExtra", "extrafont","devtools")
lapply(lib, library, character.only=TRUE)
library(readr)
library(qqman)

### read in AF_vapeR data 
all_chr_res <- lapply(1:40,function(i) readRDS(paste0("mafs/NE","chr",i,"_AFV.rds")))

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
all_plots

# Show the plots for eigenvalue 1
all_plots[[1]] + theme(axis.text.x=element_text(size=2))
ggsave("mafs/NE_eig1.png", width = 30, height =4, dpi =300, bg ="white")

rownames<-as.data.frame(rownames(pvals))
names(rownames)=c("ID")

p.df<-as.data.frame(pvals)

new_df <- data.frame(position = rownames, p.df)

split_data <- strsplit(new_df$ID, "[:-]")
names(split_data)=c("chr", "start", "end")

# Create new columns for "chromosome," "start," and "end"
new_df$chr <- sapply(split_data, function(x) x[1])
new_df$start <- as.numeric(sapply(split_data, function(x) as.numeric(x[2])))
new_df$end <- as.numeric(sapply(split_data, function(x) as.numeric(x[3])))
new_df$midpoint <- floor((new_df$start + new_df$end) / 2)


NE_AF_df<-as.data.frame(cbind(new_df$chr, new_df$midpoint, new_df$Eigenvalue_1, new_df$Eigenvalue_2, new_df$Eigenvalue_3))
names(NE_AF_df)=c("chr", "mid", "eig1", "eig2", "eig3")
head(NE_AF_df)

contigs<-read.table("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/contigs_names.txt", header=F, sep="\t")
names(contigs) = c("name1", "name2")
name1<-as.vector(contigs$name1)
name2<-as.vector(contigs$name2)

merged_df <- merge(NE_AF_df, contigs, by.x = "chr", by.y = "name1", all.x = TRUE)
NE_AF_df$chr <- merged_df$name2

write_tsv(NE_AF_df, "mafs/NE_eigenvals.tsv")


### all of the above is file prep, can run from here now the eigenvectors p-val df has been made ###

NE_AF_df<-as.data.frame(read_tsv("mafs/NE_eigenvals.tsv"))
### add SNP column for plotting happiness
length(NE_AF_df$eig3)
NE_AF_df$SNP <- 1:nrow(NE_AF_df)

manhattan(NE_AF_df, chr="chr", bp="mid", snp="SNP", p="eig1" )

library(ggrepel)

### create a list of snps in the window of interest
snpsOfInterest<-filter(don, chr=="25" & mid>24347970 & mid<25371118)

# Prepare the dataset
NE_AF_df<-as.data.frame(NE_AF_df)
don <- NE_AF_df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(mid)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(NE_AF_df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, mid) %>%
  mutate( BPcum=mid+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
  mutate( is_annotate=ifelse(-log10(eig1)>3.5, "yes", "no")) 

# Prepare X axis
axisdf <- don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

NE_AF_manhat<-ggplot(don, aes(x=BPcum, y=-log10(eig1))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  ylim(0,5) +
  
  # Add highlighted points
  geom_point(data=snpsOfInterest, color="orange", size=2) +
  
  # Custom the theme:
  theme_bw() +
  labs(y="AF_vapeR -log10(eig1)", x="") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
NE_AF_manhat
ggsave("NE_AF_manhat.png", NE_AF_manhat,width=20,height=5,dpi=600)

### okay now for the picmin plot ### 
picmin_res<-read_tsv("./Picmin/picmin_res.tsv")

### rename the chromes 

contigs<-read.table("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/contigs_names.txt", header=F, sep="\t")
names(contigs) = c("name1", "name2")
name1<-as.vector(contigs$name1)
name2<-as.vector(contigs$name2)

merged_df <- merge(picmin_res, contigs, by.x = "chr", by.y = "name1", all.x = TRUE)
picmin_res$chr <- merged_df$name2
picmin_res$chr<-as.numeric(picmin_res$chr)
picmin_res<-as.data.frame(picmin_res)
manhattan(picmin_res, chr="chr", bp="mid", snp="mid", p="picmin_p" )

## add snp column for plotting happiness
picmin_res$SNP <- 1:nrow(picmin_res)

library(ggrepel)
library(dplyr)

### create a list of snps in the window of interest
snpsOfInterest<-filter(don, chr=="25" & mid>23475000 & mid<25675000)

# Prepare the dataset
don <- picmin_res %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(mid)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(picmin_res, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, mid) %>%
  mutate( BPcum=mid+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
  mutate( is_annotate=ifelse(-log10(picmin_p)>3.5, "yes", "no")) 



# Prepare X axis
axisdf <- don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

picmin_manhat<-ggplot(don, aes(x=BPcum, y=-log10(picmin_p))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  ylim(0,5) +
  
  # Add highlighted points
  geom_point(data=snpsOfInterest, color="orange", size=2) +
  
  # Custom the theme:
  theme_bw() +
  labs(y="PicMin -log10(p)", x="") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
picmin_manhat
ggsave("picmin_manhat.png", picmin_manhat,width=20,height=5,dpi=600)

#### PBS plotting #### 

TYD_pbs<-read.table("TYD_chr.pbs.fst.txt", sep = '\t', header=F)
TYE_pbs<-read.table("TYE_chr.pbs.fst.txt", sep = '\t', header=F)
WEA_pbs<-read.table("WEA_chr.pbs.fst.txt", sep = '\t', header=F)

names(TYD_pbs) = c("region",  "chr", "midPos",  "Nsites",  "Fst01",   "Fst02",   "Fst12", "PBS0", "PBS1", "PBS2")
TYD_pbs$pval<-qvalue::empPvals(stat = TYD_pbs$PBS0,stat0 = TYD_pbs$PBS0)

names(TYE_pbs) = c("region",  "chr", "midPos",  "Nsites",  "Fst01",   "Fst02",   "Fst12", "PBS0", "PBS1", "PBS2")
TYE_pbs$pval<-qvalue::empPvals(stat = TYE_pbs$PBS0,stat0 = TYE_pbs$PBS0)

names(WEA_pbs) = c("region",  "chr", "midPos",  "Nsites",  "Fst01",   "Fst02",   "Fst12", "PBS0", "PBS1", "PBS2")
WEA_pbs$pval<-qvalue::empPvals(stat = WEA_pbs$PBS0,stat0 = WEA_pbs$PBS0)

### TYD ### 

###rename chrs ### 

contigs<-read.table("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/contigs_names.txt", header=F, sep="\t")
names(contigs) = c("name1", "name2")
name1<-as.vector(contigs$name1)
name2<-as.vector(contigs$name2)

merged_df <- merge(TYD_pbs, contigs, by.x = "chr", by.y = "name1", all.x = TRUE)
TYD_pbs$chr <- merged_df$name2
TYD_pbs$chr<-as.numeric(TYD_pbs$chr)
TYD_pbs<-as.data.frame(TYD_pbs)
TYD_pbs$chr<-as.numeric(TYD_pbs$chr)

### create a list of snps in the window of interest
snpsOfInterest<-filter(TYD_pbs, chr=="25" & midPos>24425000 & midPos<25675000) ### taken from pbs windows where p<0.05 in TYD
manhattan(TYD_pbs, chr="chr", bp="midPos", snp="midPos", p="pval", highlight = snpsOfInterest )

# ## add snp column for plotting happiness
# TYD_pbs$SNP <- 1:nrow(TYD_pbs)
# TYD_pbs<-as.data.frame(TYD_pbs)
# 
# library(ggrepel)
# library(dplyr)
# 
# 
# # Prepare the dataset
# don <- TYD_pbs %>% 
#   
#   # Compute chromosome size
#   group_by(chr) %>% 
#   summarise(chr_len=max(midPos)) %>% 
#   
#   # make chr_len numeric
#   don$chr_len=as.numeric(don$chr_len) %>%
# 
#   # Calculate cumulative position of each chromosome
#   group_by(chr) %>%
#   mutate(tot=cumsum(chr_len)-chr_len) %>%
#   ungroup() %>%
#   dplyr::select(-chr_len)# %>%
#   
#   # Add this info to the initial dataset
#   left_join(TYD_pbs, ., by=c("chr"="chr")) %>%
#   
#   # Add a cumulative position of each SNP
#   arrange(chr, midPos) %>%
#   mutate( BPcum=midPos+tot) %>%
#   
#   # Add highlight and annotation information
#   mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
#   mutate( is_annotate=ifelse(-log10(pval)>3.5, "yes", "no")) 
# 
# 
# # Prepare X axis
# axisdf <- don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
# 
# TYD_PBS_manhat<-ggplot(don, aes(x=BPcum, y=-log10(pval))) +
#   
#   # Show all points
#   geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
#   scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
#   
#   # custom X axis:
# scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
# scale_y_continuous(expand = c(0, 0) )# +  # remove space between plot area and x axis
# #  ylim(0,5) +
#   
#   # Add highlighted points
#   geom_point(data=snpsOfInterest, color="orange", size=2) +
#   
#   # Custom the theme:
# #  theme_bw() +
# #  labs(y="TYD PBS -log10(p)", x="") +
# #  theme( 
# #    legend.position="none",
# #    panel.border = element_blank(),
# #    panel.grid.major.x = element_blank(),
# #    panel.grid.minor.x = element_blank()
#   )
# TYD_PBS_manhat
# ggsave("picmin_manhat.png", picmin_manhat,width=20,height=5,dpi=600)





### PCAngsd #### 

# library(RcppCNPy)
# source("./pcangsd_plotting/individual_pca_functions.R")
# 
# C <- as.matrix(read.table("pcangsd/pcangsd_NE_snp_selection.cov")) # Reads in estimated covariance matrix
# D <- npyLoad("pcangsd/pcangsd_NE_snp_selection.selection.npy") # Reads PC based selection statistics
# 
# samples <- read.csv("./pcangsd/sample_list_wales.csv", header = F)
# names(samples) =c ("ind", "pop", "path")
# ind <- as.matrix(samples$ind)
# pop <- as.matrix(samples$pop)
# 
# # Obtain p-values from selection scan
# pval <- pchisq(D, 1, lower.tail=FALSE)
# 
# allsites<-read.table("./pcangsd/a.pos", colC=c("factor","integer"),sep="\t")
# used<-read.table("./pcangsd/pcangsd_NE_snp_selection.sites",colC=c("integer"),sep="\t")
# names(used)<-c("inc")
# #sites<-read.table("pcangsd_wales_selection.sites",colC=c("factor","integer"),sep="\t")
# names(allsites)<-c("chr","pos")
# J<-cbind.data.frame(allsites$chr, allsites$pos, used$inc)
# K<- data.frame(subset(J, used$inc =="1"))
# names(K) = c("chr", "pos", "inc")
# sites<-cbind.data.frame(K$chr, K$pos)
# names(sites) = c("chr", "pos")
# 
# # now combine the pvalues calculated from the test statistics with the filtered positions
# #each column of the pval is a PC axis (e.g. pval[,1]=PC1)
# 
# pval<-as.data.frame(pval)
#pca_matrix <- cbind.data.frame(sites$chr, sites$pos, pval$V1, pval$V2, pval[,3])# pval[,4], pval[,5], pval[,6])
#names(pca_matrix) <- c("chr", "pos", "pc1_pval", "pc2_pval", "pc3_pval")#, "pc4_pval", "pc5_pval", "pc6_pval")
# pca_matrix.df <- as.data.frame(pca_matrix)
# 
# contigs<-read.table("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/contigs_names.txt", header=F, sep="\t")
# names(contigs) = c("name1", "name2")
# name1<-as.vector(contigs$name1)
# name2<-as.vector(contigs$name2)
# 
# merged_df <- merge(pca_matrix.df, contigs, by.x = "chr", by.y = "name1", all.x = TRUE)
# pca_matrix.df$chr <- merged_df$name2

#write_csv(pca_matrix.df, "./pcangsd/NE_full_pca_matrix.df.csv")



### above is just dataprocessing ###
### read back in df ###
pca_matrix.df<-read_csv("./pcangsd/NE_full_pca_matrix.df.csv")
head(pca_matrix.df)

pca_matrix.df$SNP <- 1:nrow(pca_matrix.df)
head(pca_matrix.df)

pca_matrix.df <- pca_matrix.df %>%
  filter(-log10(pc3_pval)>1)


#manhattan(pca_matrix.df, chr="chr", bp="pos", snp="SNP", p="pc3_pval" )

library(ggrepel)

### create a list of snps in the window of interest
snpsOfInterest<-filter(don, chr=="25" & pos>24347970 & pos<25371118)

# Prepare the dataset
pca_matrix.df<-as.data.frame(pca_matrix.df)
don <- pca_matrix.df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(pca_matrix.df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( BPcum=pos+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
  mutate( is_annotate=ifelse(-log10(pc3_pval)>3.5, "yes", "no")) 

# Prepare X axis
axisdf <- don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

NE_pcangsd_manhat<-ggplot(don, aes(x=BPcum, y=-log10(pc3_pval))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "palegreen"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  ylim(0,8) +
  
  # Add highlighted points
  geom_point(data=snpsOfInterest, color="orange", size=2) +
  
  # Custom the theme:
  theme_bw() +
  labs(y="NE pcangsd -log10(eig1)", x="") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
NE_pcangsd_manhat
ggsave("pcangsd_manhat.png", NE_pcangsd_manhat,width=20,height=5,dpi=600)

library(cowplot)
com_plot<-plot_grid(picmin_manhat, NE_AF_manhat, NE_pcangsd_manhat,
          ncol=1)
ggsave("combined_plot.png", com_plot, width=10, height = 8)

