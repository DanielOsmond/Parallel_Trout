setwd("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/sweepfinder/output")
library(readr)
library(tidyverse)

TYD<-as.data.frame(read.table("TYD_25_filt.sw.out", header= T, sep = "\t"))
TYA<-as.data.frame(read.table("TYA_25_filt.sw.out", header= T, sep = "\t"))
WEF<-as.data.frame(read.table("WEF_25_filt.sw.out", header= T, sep = "\t"))
WEA<-as.data.frame(read.table("WEA_25_filt.sw.out", header= T, sep = "\t"))
TYE<-as.data.frame(read.table("TYE_25_filt.sw.out", header= T, sep = "\t"))
head(TYD)

ggplot(TYA) +
  geom_point(aes(x = location, y = LR)) +
  theme_bw()


TYD_filt<-filter(TYD, alpha < 0.05)

TYE_plot <-ggplot(TYE) +
  geom_point(aes(x = location, y = LR)) +
  theme_bw()

TYA_plot<-ggplot(TYA) +
  geom_point(aes(x = location, y = LR)) +
  theme_bw()

TYD_plot<-ggplot(TYD) +
  geom_point(aes(x = location, y = LR)) +
  theme_bw()

WEA_plot<-ggplot(WEA) +
  geom_point(aes(x = location, y = LR)) +
  theme_bw()

WEF_plot<-ggplot(WEF) +
  geom_point(aes(x = location, y = LR)) +
  theme_bw()

ne_sweep_25<-plot_grid(TYA_plot, TYD_plot, TYE_plot, WEA_plot, WEF_plot,
          ncol = 1)

ggsave("ne_sweep_25.png", width =6, height =10)

max(TYA$LR)



### looped code from openai ### 
# Initialize an empty list to store data frames
data_frames <- list()

# Loop through the files
for (i in 1:40) {
  # Construct the file name
  file_name <- paste("HAY1_", i, "_filt.sw.out", sep = "")
  
  # Read the data into a data frame
  data <- as.data.frame(read.table(file_name, header = TRUE, sep = "\t"))
  
  # Store the data frame in the list
  data_frames[[i]] <- data
  
  # Create a plot with the appropriate name
  plot_name <- paste("HAY_", i, "_plot", sep = "")
  
  # Create the plot
  assign(plot_name, ggplot(data) +
           geom_point(aes(x = location, y = LR)) +
           theme_bw()
  )
  
  # Save the plot as a file (optional)
  ggsave(paste(plot_name, ".png", sep = ""), plot = get(plot_name), width = 8, height = 6)
}

#########################################################
#### merged SW outputs ##################################
#########################################################

library(readr)
library(ggplot2)
library(dplyr)
library(data.table)

## AGA ##

# AGA<-as.data.frame(read.table("../AGA_merged.txt", header=T, sep = "\t"))
# names(AGA)=c("pos", "CLR", "alpha", "chr", "pop")
# AGA$CLR<-as.numeric(AGA$CLR)
# AGA$pos<-as.numeric(AGA$pos)
# summary(AGA$CLR)
# q99.9<-quantile(AGA$CLR, probs=c(0.999), na.rm =T)
# 
# # Prepare the dataset
# AGA.df<-na.omit(as.data.frame(AGA))
# don <- AGA.df %>% 
#   
#   # Compute chromosome size
#   group_by(chr) %>% 
#   summarise(chr_len=max(pos)) %>% 
#   
#   # Calculate cumulative position of each chromosome
#   mutate(tot=cumsum(chr_len)-chr_len) %>%
#   dplyr::select(-chr_len) %>%
#   
#   # Add this info to the initial dataset
#   left_join(AGA.df, ., by=c("chr"="chr")) %>%
#   
#   # Add a cumulative position of each SNP
#   arrange(chr, pos) %>%
#   mutate( BPcum=pos+tot) #%>%
#   
#   # Add highlight and annotation information
# #  mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
# #  mutate( is_annotate=ifelse(-log10(pc3_pval)>3.5, "yes", "no")) 
# 
# # Prepare X axis
# axisdf <- don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
# 
# AGA_sweep_manhat<-ggplot(don, aes(x=BPcum, y= CLR)) +
#   
#   # Show all points
#   geom_point( aes(color=as.factor(chr)), alpha=0.8, size=0.3) +
#   scale_color_manual(values = rep(c("grey", "palegreen"), 22 )) +
#   
#   # custom X axis:
#   scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
#   #scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
#  # ylim(0,8) +
#   
#   # Add highlighted points
# #  geom_point(data=snpsOfInterest, color="orange", size=2) +
#   
#   # Custom the theme:
#   theme_bw() +
# #  labs(y="NE pcangsd -log10(eig1)", x="") +
#   theme( 
#     legend.position="none",
#     panel.border = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank()
#   )
# AGA_sweep_manhat
# ggsave("AGA_sweep.png", AGA_sweep_manhat,width=20,height=5,dpi=600)
# 
# 
# ## Hayle ##
# 
# HAY<-as.data.frame(read.table("../HAY1_merged.txt", header=T, sep = "\t"))
# names(HAY)=c("pos", "CLR", "alpha", "chr", "pop")
# HAY$CLR<-as.numeric(HAY$CLR)
# HAY$pos<-as.numeric(HAY$pos)
# 
# q99.9<-quantile(HAY$CLR, probs=c(0.999), na.rm =T)
# 
# # Prepare the dataset
# HAY.df<-na.omit(as.data.frame(HAY))
# hon <- HAY.df %>% 
#   
#   # Compute chromosome size
#   group_by(chr) %>% 
#   summarise(chr_len=max(pos)) %>% 
#   
#   # Calculate cumulative position of each chromosome
#   mutate(tot=cumsum(chr_len)-chr_len) %>%
#   dplyr::select(-chr_len) %>%
#   
#   # Add this info to the initial dataset
#   left_join(HAY.df, ., by=c("chr"="chr")) %>%
#   
#   # Add a cumulative position of each SNP
#   arrange(chr, pos) %>%
#   mutate( BPcum=pos+tot) #%>%
# 
# # Add highlight and annotation information
# #  mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
# #  mutate( is_annotate=ifelse(-log10(pc3_pval)>3.5, "yes", "no")) 
# 
# # Prepare X axis
# axisdf <- hon %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
# 
# HAY_sweep_manhat<-ggplot(hon, aes(x=BPcum, y= CLR)) +
#   
#   # Show all points
#   geom_point( aes(color=as.factor(chr)), alpha=0.8, size=0.3) +
#   scale_color_manual(values = rep(c("grey", "palegreen"), 22 )) +
#   
#   # custom X axis:
#   scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
#   #scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
#   # ylim(0,8) +
#   
#   # Add highlighted points
#   #  geom_point(data=snpsOfInterest, color="orange", size=2) +
#   
#   # Custom the theme:
#   theme_bw() +
#   #  labs(y="NE pcangsd -log10(eig1)", x="") +
#   theme( 
#     legend.position="none",
#     panel.border = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank()
#   )
# HAY_sweep_manhat
# ggsave("HAY_sweep.png", HAY_sweep_manhat,width=20,height=5,dpi=600)
# 
subsetDataFrameByPercentile <- function(data, column, percentile) {
  q <- quantile(data[[column]], probs = percentile, na.rm = TRUE)
  subset_data <- na.omit(data[data[[column]] >= q, ])
  return(subset_data)
}
# 
# outliers_HAY <- subsetDataFrameByPercentile(HAY, "CLR", 0.995)
# outliers_AGA <- subsetDataFrameByPercentile(AGA, "CLR", 0.995)

### let's pump those outliers ###

AGA<-as.data.frame(read.table("../AGA_merged.txt", header=T, sep = "\t"))
names(AGA)=c("pos", "CLR", "alpha", "chr", "pop")
AGA$CLR<-as.numeric(AGA$CLR)
AGA$pos<-as.numeric(AGA$pos)
AGA$ID<-paste(AGA$chr, ":", AGA$pos)
outliers_AGA <- subsetDataFrameByPercentile(AGA, "CLR", 0.995)

AGB<-as.data.frame(read.table("../AGB_merged.txt", header=T, sep = "\t"))
names(AGB)=c("pos", "CLR", "alpha", "chr", "pop")
AGB$CLR<-as.numeric(AGB$CLR)
AGB$pos<-as.numeric(AGB$pos)
AGB$ID<-paste(AGB$chr, ":", AGB$pos)
outliers_AGB <- subsetDataFrameByPercentile(AGB, "CLR", 0.995)

AVB<-as.data.frame(read.table("../AVB_merged.txt", header=T, sep = "\t"))
names(AVB)=c("pos", "CLR", "alpha", "chr", "pop")
AVB$CLR<-as.numeric(AVB$CLR)
AVB$pos<-as.numeric(AVB$pos)
AVB$ID<-paste(AVB$chr, ":", AVB$pos)
outliers_AVB <- subsetDataFrameByPercentile(AVB, "CLR", 0.995)

# AVD<-as.data.frame(read.table("../AVD_merged.txt", header=T, sep = "\t"))
# names(AVD)=c("pos", "CLR", "alpha", "chr", "pop")
# AVD$CLR<-as.numeric(AVD$CLR)
# AVD$pos<-as.numeric(AVD$pos)
# AVD$ID<-paste(AVD$chr, ":", AVD$pos)
# outliers_AGB <- subsetDataFrameByPercentile(AVD, "CLR", 0.995)

# CAM<-as.data.frame(read.table("../CAM_merged.txt", header=T, sep = "\t"))
# names(CAM)=c("pos", "CLR", "alpha", "chr", "pop")
# CAM$CLR<-as.numeric(CAM$CLR)
# CAM$pos<-as.numeric(CAM$pos)
# CAM$ID<-paste(CAM$chr, ":", CAM$pos)
# outliers_CAM <- subsetDataFrameByPercentile(CAM, "CLR", 0.995)

CWA<-as.data.frame(read.table("../CWA_merged.txt", header=T, sep = "\t"))
names(CWA)=c("pos", "CLR", "alpha", "chr", "pop")
CWA$CLR<-as.numeric(CWA$CLR)
CWA$pos<-as.numeric(CWA$pos)
CWA$ID<-paste(CWA$chr, ":", CWA$pos)
outliers_CWA <- subsetDataFrameByPercentile(CWA, "CLR", 0.995)

CWB<-as.data.frame(read.table("../CWB_merged.txt", header=T, sep = "\t"))
names(CWB)=c("pos", "CLR", "alpha", "chr", "pop")
CWB$CLR<-as.numeric(CWB$CLR)
CWB$pos<-as.numeric(CWB$pos)
CWB$ID<-paste(CWB$chr, ":", CWB$pos)
outliers_CWB <- subsetDataFrameByPercentile(CWB, "CLR", 0.995)

HAY<-as.data.frame(read.table("../HAY1_merged.txt", header=T, sep = "\t"))
names(HAY)=c("pos", "CLR", "alpha", "chr", "pop")
HAY$CLR<-as.numeric(HAY$CLR)
HAY$pos<-as.numeric(HAY$pos)
HAY$ID<-paste(HAY$chr, ":", HAY$pos)
outliers_HAY <- subsetDataFrameByPercentile(HAY, "CLR", 0.995)

# NEV<-as.data.frame(read.table("../NEV2_merged.txt", header=T, sep = "\t"))
# names(NEV)=c("pos", "CLR", "alpha", "chr", "pop")
# NEV$CLR<-as.numeric(NEV$CLR)
# NEV$pos<-as.numeric(NEV$pos)
# NEV$ID<-paste(NEV$chr, ":", NEV$pos)
# outliers_NEV <- subsetDataFrameByPercentile(NEV, "CLR", 0.995)

# OWE<-as.data.frame(read.table("../OWE_merged.txt", header=T, sep = "\t"))
# names(OWE)=c("pos", "CLR", "alpha", "chr", "pop")
# OWE$CLR<-as.numeric(OWE$CLR)
# OWE$pos<-as.numeric(OWE$pos)
# OWE$ID<-paste(OWE$chr, ":", OWE$pos)
# outliers_OWE <- subsetDataFrameByPercentile(OWE, "CLR", 0.995)

TYA<-as.data.frame(read.table("../TYA_merged.txt", header=T, sep = "\t"))
names(TYA)=c("pos", "CLR", "alpha", "chr", "pop")
TYA$CLR<-as.numeric(TYA$CLR)
TYA$pos<-as.numeric(TYA$pos)
TYA$ID<-paste(TYA$chr, ":", TYA$pos)
outliers_TYA <- subsetDataFrameByPercentile(TYA, "CLR", 0.995)

TYD<-as.data.frame(read.table("../TYD_merged.txt", header=T, sep = "\t"))
names(TYD)=c("pos", "CLR", "alpha", "chr", "pop")
TYD$CLR<-as.numeric(TYD$CLR)
TYD$pos<-as.numeric(TYD$pos)
TYD$ID<-paste(TYD$chr, ":", TYD$pos)
outliers_TYD <- subsetDataFrameByPercentile(TYD, "CLR", 0.995)

TYE<-as.data.frame(read.table("../TYE_merged.txt", header=T, sep = "\t"))
names(TYE)=c("pos", "CLR", "alpha", "chr", "pop")
TYE$CLR<-as.numeric(TYE$CLR)
TYE$pos<-as.numeric(TYE$pos)
TYE$ID<-paste(TYE$chr, ":", TYE$pos)
outliers_TYE <- subsetDataFrameByPercentile(TYE, "CLR", 0.995)

WEA<-as.data.frame(read.table("../WEA_merged.txt", header=T, sep = "\t"))
names(WEA)=c("pos", "CLR", "alpha", "chr", "pop")
WEA$CLR<-as.numeric(WEA$CLR)
WEA$pos<-as.numeric(WEA$pos)
WEA$ID<-paste(WEA$chr, ":", WEA$pos)
outliers_WEA <- subsetDataFrameByPercentile(WEA, "CLR", 0.995)

WEF<-as.data.frame(read.table("../WEF_merged.txt", header=T, sep = "\t"))
names(WEF)=c("pos", "CLR", "alpha", "chr", "pop")
WEF$CLR<-as.numeric(WEF$CLR)
WEF$pos<-as.numeric(WEF$pos)
WEF$ID<-paste(WEF$chr, ":", WEF$pos)
outliers_WEF <- subsetDataFrameByPercentile(WEF, "CLR", 0.995)

YHA<-as.data.frame(read.table("../YHA_merged.txt", header=T, sep = "\t"))
names(YHA)=c("pos", "CLR", "alpha", "chr", "pop")
YHA$CLR<-as.numeric(YHA$CLR)
YHA$pos<-as.numeric(YHA$pos)
YHA$ID<-paste(YHA$chr, ":", YHA$pos)
outliers_YHA <- subsetDataFrameByPercentile(YHA, "CLR", 0.995)

YHE<-as.data.frame(read.table("../YHE_merged.txt", header=T, sep = "\t"))
names(YHE)=c("pos", "CLR", "alpha", "chr", "pop")
YHE$CLR<-as.numeric(YHE$CLR)
YHE$pos<-as.numeric(YHE$pos)
YHE$ID<-paste(YHE$chr, ":", YHE$pos)
outliers_YHE <- subsetDataFrameByPercentile(YHE, "CLR", 0.995)

#merged_sw<-setDT(AGA)[ID %chin%AGB$ID]
### this approach cuts about 2/3 of the sites, I guess because the positions don't match exactly ### 
### I think best approach is to overlap each PBS stat with each SW outlier to prevent 

library(GenomicRanges)

### read in my PBS matrix ### 
pbs<-read_tsv("../../sfs/pbs/gene_annotated_pbs.txt")
pbs<-subset(pbs, AGA_pval<0.01)

pbs_gr<-makeGRangesFromDataFrame(pbs,
                              keep.extra.columns = TRUE,
                              seqinfo = NULL,
                              seqnames.field=c("chromosome_name"),
                              start.field = "start", 
                              end.field = "end",
                              ignore.strand = TRUE,
                              na.rm=TRUE)

AGA_gr<-makeGRangesFromDataFrame(outliers_AGA,
                              keep.extra.columns = TRUE,
                              seqinfo = NULL,
                              seqnames.field=c("chr"),
                              start.field = "pos", 
                              end.field = "pos",
                              ignore.strand = TRUE,
                              na.rm=TRUE)

AGA_sub<-subsetByOverlaps(pbs_gr, AGA_gr)
AGA_sub2<-subsetByOverlaps(AGA_gr, pbs_gr)
AGA_overlap<-as.data.frame(AGA_sub)
AGA_overlap2<-as.data.frame(AGA_sub2)
write_tsv(AGA_overlap, "AGA_pbs_sw_overlap_genes.txt")
write_tsv(AGA_overlap2, "AGA_sw_pbs_overlap_genes.txt")
AGA_overlap_pbs<-as.data.frame(cbind(AGA_overlap$seqnames, AGA_overlap$midPos, AGA_overlap$AGA_pval, AGA_overlap$ensembl_gene_id, AGA_overlap$description, AGA_overlap$external_gene_name))
names(AGA_overlap_pbs)=c("chr", "midpos", "emppval", "ensembl_geneID", "description", "gene_name")
AGA_overlap_pbs<-subset(AGA_overlap_pbs, emppval<0.01)


# Prepare the dataset
AGA.df<-na.omit(as.data.frame(AGA))
don_AGA <- AGA.df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(AGA.df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( BPcum=pos+tot) #%>%

# Add highlight and annotation information
#  mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
#  mutate( is_annotate=ifelse(-log10(pc3_pval)>3.5, "yes", "no")) 

# Prepare X axis
axisdf <- don_AGA %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

highlight_AGA <- inner_join(don_AGA, AGA_overlap2, by = "ID")
names(highlight_AGA)=c("pos", "CLR", "alpha", "chr", "pop", "ID", "tot", "BPcum", "seqnames", "start", "end", "width", "strand", "CLR.y", "alpha.y", "pop.y")

AGA_sweep_manhat<-ggplot(don_AGA, aes(x=BPcum, y= CLR)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=0.3) +
  scale_color_manual(values = rep(c("grey", "palegreen"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  # ylim(0,8) +
  
  # Add highlighted points
 geom_point(data=highlight_AGA, color="red", size=0.5) +
  
  # Custom the theme:
  theme_bw() +
   labs(y="CLR", x="AGA") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
#AGA_sweep_manhat
#ggsave("AGA_sweep_pbs.png", AGA_sweep_manhat,width=20,height=5,dpi=600)

## CWB ####################################################################################################################
#############################################################################################################################
###########################################################################################################################

pbs<-read_tsv("../../sfs/pbs/gene_annotated_pbs.txt")
pbs<-subset(pbs, CWB_pval<0.01)

pbs_gr<-makeGRangesFromDataFrame(pbs,
                                 keep.extra.columns = TRUE,
                                 seqinfo = NULL,
                                 seqnames.field=c("chromosome_name"),
                                 start.field = "start", 
                                 end.field = "end",
                                 ignore.strand = TRUE,
                                 na.rm=TRUE)

CWB_gr<-makeGRangesFromDataFrame(outliers_CWB,
                                 keep.extra.columns = TRUE,
                                 seqinfo = NULL,
                                 seqnames.field=c("chr"),
                                 start.field = "pos", 
                                 end.field = "pos",
                                 ignore.strand = TRUE,
                                 na.rm=TRUE)

CWB_sub<-subsetByOverlaps(pbs_gr, CWB_gr)
CWB_sub2<-subsetByOverlaps(CWB_gr, pbs_gr)
CWB_overlap<-as.data.frame(CWB_sub)
CWB_overlap2<-as.data.frame(CWB_sub2)
write_tsv(CWB_overlap, "CWB_pbs_sw_overlap_genes.txt")
write_tsv(CWB_overlap2, "CWB_sw_pbs_overlap_genes.txt")
CWB_overlap_pbs<-as.data.frame(cbind(CWB_overlap$seqnames, CWB_overlap$midPos, CWB_overlap$CWB_pval, CWB_overlap$ensembl_gene_id, CWB_overlap$description, CWB_overlap$external_gene_name))
names(CWB_overlap_pbs)=c("chr", "midpos", "emppval", "ensembl_geneID", "description", "gene_name")
CWB_overlap_pbs<-subset(CWB_overlap_pbs, emppval<0.01)

# Prepare the dataset
CWB.df<-na.omit(as.data.frame(CWB))
don_CWB <- CWB.df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(CWB.df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( BPcum=pos+tot) #%>%

# Add highlight_CWB and annotation information
#  mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
#  mutate( is_annotate=ifelse(-log10(pc3_pval)>3.5, "yes", "no")) 

# Prepare X axis
axisdf <- don_CWB %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

highlight_CWB <- inner_join(don_CWB, CWB_overlap2, by = "ID")
names(highlight_CWB)=c("pos", "CLR", "alpha", "chr", "pop", "ID", "tot", "BPcum", "seqnames", "start", "end", "width", "strand", "CLR.y", "alpha.y", "pop.y")

CWB_sweep_manhat<-ggplot(don_CWB, aes(x=BPcum, y= CLR)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=0.3) +
  scale_color_manual(values = rep(c("grey", "palegreen"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  # ylim(0,8) +
  
  # Add highlighted points
  geom_point(data=highlight_CWB, color="red", size=0.3) +
  
  # Custom the theme:
  theme_bw() +
   labs(y="CLR", x="CWB") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
#CWB_sweep_manhat
#ggsave("CWB_sweep_pbs.png", CWB_sweep_manhat,width=20,height=5,dpi=600)

#################################################################################
###################### HAYLE ####################################################
#################################################################################

pbs<-read_tsv("../../sfs/pbs/gene_annotated_pbs.txt")
pbs<-subset(pbs, HAY_pval<0.01)

pbs_gr<-makeGRangesFromDataFrame(pbs,
                                 keep.extra.columns = TRUE,
                                 seqinfo = NULL,
                                 seqnames.field=c("chromosome_name"),
                                 start.field = "start", 
                                 end.field = "end",
                                 ignore.strand = TRUE,
                                 na.rm=TRUE)

HAY_gr<-makeGRangesFromDataFrame(outliers_HAY,
                                 keep.extra.columns = TRUE,
                                 seqinfo = NULL,
                                 seqnames.field=c("chr"),
                                 start.field = "pos", 
                                 end.field = "pos",
                                 ignore.strand = TRUE,
                                 na.rm=TRUE)

HAY_sub<-subsetByOverlaps(pbs_gr, HAY_gr)
HAY_sub2<-subsetByOverlaps(HAY_gr, pbs_gr)
HAY_overlap<-as.data.frame(HAY_sub)
HAY_overlap2<-as.data.frame(HAY_sub2)
write_tsv(HAY_overlap, "HAY_pbs_sw_overlap_genes.txt")
write_tsv(HAY_overlap2, "HAY_sw_pbs_overlap_genes.txt")
HAY_overlap_pbs<-as.data.frame(cbind(HAY_overlap$seqnames, HAY_overlap$midPos, HAY_overlap$HAY_pval, HAY_overlap$ensembl_gene_id, HAY_overlap$description, HAY_overlap$external_gene_name))
names(HAY_overlap_pbs)=c("chr", "midpos", "emppval", "ensembl_geneID", "description", "gene_name")
HAY_overlap_pbs<-subset(HAY_overlap_pbs, emppval<0.01)

# Prepare the dataset
HAY.df<-na.omit(as.data.frame(HAY))
don_HAY <- HAY.df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(HAY.df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( BPcum=pos+tot) #%>%

# Add highlight_HAY and annotation information
#  mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
#  mutate( is_annotate=ifelse(-log10(pc3_pval)>3.5, "yes", "no")) 

# Prepare X axis
axisdf <- don_HAY %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

highlight_HAY <- inner_join(don_HAY, HAY_overlap2, by = "ID")
names(highlight_HAY)=c("pos", "CLR", "alpha", "chr", "pop", "ID", "tot", "BPcum", "seqnames", "start", "end", "width", "strand", "CLR.y", "alpha.y", "pop.y")

HAY_sweep_manhat<-ggplot(don_HAY, aes(x=BPcum, y= CLR)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=0.3) +
  scale_color_manual(values = rep(c("grey", "palegreen"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  # ylim(0,8) +
  
  # Add highlighted points
  geom_point(data=highlight_HAY, color="red", size=0.3) +
  
  # Custom the theme:
  theme_bw() +
  labs(y="CLR", x="HAY") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
#HAY_sweep_manhat
#ggsave("HAY_sweep_pbs.png", HAY_sweep_manhat,width=20,height=5,dpi=600)

#################################################################################
###################### Avoca ####################################################
#################################################################################

pbs<-read_tsv("../../sfs/pbs/gene_annotated_pbs.txt")
pbs<-subset(pbs, AVB_pval<0.01)

pbs_gr<-makeGRangesFromDataFrame(pbs,
                                 keep.extra.columns = TRUE,
                                 seqinfo = NULL,
                                 seqnames.field=c("chromosome_name"),
                                 start.field = "start", 
                                 end.field = "end",
                                 ignore.strand = TRUE,
                                 na.rm=TRUE)

AVB_gr<-makeGRangesFromDataFrame(outliers_AVB,
                                 keep.extra.columns = TRUE,
                                 seqinfo = NULL,
                                 seqnames.field=c("chr"),
                                 start.field = "pos", 
                                 end.field = "pos",
                                 ignore.strand = TRUE,
                                 na.rm=TRUE)

AVB_sub<-subsetByOverlaps(pbs_gr, AVB_gr)
AVB_sub2<-subsetByOverlaps(AVB_gr, pbs_gr)
AVB_overlap<-as.data.frame(AVB_sub)
AVB_overlap2<-as.data.frame(AVB_sub2)
write_tsv(AVB_overlap, "AVB_pbs_sw_overlap_genes.txt")
write_tsv(AVB_overlap2, "AVB_sw_pbs_overlap_genes.txt")
AVB_overlap_pbs<-as.data.frame(cbind(AVB_overlap$seqnames, AVB_overlap$midPos, AVB_overlap$AVB_pval, AVB_overlap$ensembl_gene_id, AVB_overlap$description, AVB_overlap$external_gene_name))
names(AVB_overlap_pbs)=c("chr", "midpos", "emppval", "ensembl_geneID", "description", "gene_name")
AVB_overlap_pbs<-subset(AVB_overlap_pbs, emppval<0.01)

# Prepare the dataset
AVB.df<-na.omit(as.data.frame(AVB))
don_AVB <- AVB.df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(AVB.df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( BPcum=pos+tot) #%>%

# Add highlight_AVB and annotation information
#  mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
#  mutate( is_annotate=ifelse(-log10(pc3_pval)>3.5, "yes", "no")) 

# Prepare X axis
axisdf <- don_AVB %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

highlight_AVB <- inner_join(don_AVB, AVB_overlap2, by = "ID")
names(highlight_AVB)=c("pos", "CLR", "alpha", "chr", "pop", "ID", "tot", "BPcum", "seqnames", "start", "end", "width", "strand", "CLR.y", "alpha.y", "pop.y")

AVB_sweep_manhat<-ggplot(don_AVB, aes(x=BPcum, y= CLR)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=0.3) +
  scale_color_manual(values = rep(c("grey", "palegreen"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  # ylim(0,8) +
  
  # Add highlight_AVBed points
  geom_point(data=highlight_AVB, color="red", size=0.3) +
  
  # Custom the theme:
  theme_bw() +
  labs(y="CLR", x="AVB") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
#AVB_sweep_manhat
#ggsave("AVB_sweep_pbs.png", AVB_sweep_manhat,width=20,height=5,dpi=600)

#################################################################################
###################### TYD ######################################################
#################################################################################

pbs<-read_tsv("../../sfs/pbs/gene_annotated_pbs.txt")
pbs<-subset(pbs, TYD_pval<0.01)

pbs_gr<-makeGRangesFromDataFrame(pbs,
                                 keep.extra.columns = TRUE,
                                 seqinfo = NULL,
                                 seqnames.field=c("chromosome_name"),
                                 start.field = "start", 
                                 end.field = "end",
                                 ignore.strand = TRUE,
                                 na.rm=TRUE)

TYD_gr<-makeGRangesFromDataFrame(outliers_TYD,
                                 keep.extra.columns = TRUE,
                                 seqinfo = NULL,
                                 seqnames.field=c("chr"),
                                 start.field = "pos", 
                                 end.field = "pos",
                                 ignore.strand = TRUE,
                                 na.rm=TRUE)

TYD_sub<-subsetByOverlaps(pbs_gr, TYD_gr)
TYD_sub2<-subsetByOverlaps(TYD_gr, pbs_gr)
TYD_overlap<-as.data.frame(TYD_sub)
TYD_overlap2<-as.data.frame(TYD_sub2)
write_tsv(TYD_overlap, "TYD_pbs_sw_overlap_genes.txt")
write_tsv(TYD_overlap2, "TYD_sw_pbs_overlap_genes.txt")
TYD_overlap_pbs<-as.data.frame(cbind(TYD_overlap$seqnames, TYD_overlap$midPos, TYD_overlap$TYD_pval, TYD_overlap$ensembl_gene_id, TYD_overlap$description, TYD_overlap$external_gene_name))
names(TYD_overlap_pbs)=c("chr", "midpos", "emppval", "ensembl_geneID", "description", "gene_name")
TYD_overlap_pbs<-subset(TYD_overlap_pbs, emppval<0.01)

# Prepare the dataset
TYD.df<-na.omit(as.data.frame(TYD))
don_TYD <- TYD.df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(TYD.df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( BPcum=pos+tot) #%>%

# Add highlight_TYD and annotation information
#  mutate( is_highlight_TYD=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
#  mutate( is_annotate=ifelse(-log10(pc3_pval)>3.5, "yes", "no")) 

# Prepare X axis
axisdf <- don_TYD %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

highlight_TYD <- inner_join(don_TYD, TYD_overlap2, by = "ID")
names(highlight_TYD)=c("pos", "CLR", "alpha", "chr", "pop", "ID", "tot", "BPcum", "seqnames", "start", "end", "width", "strand", "CLR.y", "alpha.y", "pop.y")

TYD_sweep_manhat<-ggplot(don_TYD, aes(x=BPcum, y= CLR)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=0.3) +
  scale_color_manual(values = rep(c("grey", "palegreen"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  # ylim(0,8) +
  
  # Add highlight_TYDed points
  geom_point(data=highlight_TYD, color="red", size=0.3) +
  
  # Custom the theme:
  theme_bw() +
  labs(y="CLR", x="TYD") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
#TYD_sweep_manhat
#ggsave("TYD_sweep_pbs.png", TYD_sweep_manhat,width=20,height=5,dpi=600)

#################################################################################
###################### TYE ######################################################
#################################################################################

pbs<-read_tsv("../../sfs/pbs/gene_annotated_pbs.txt")
pbs<-subset(pbs, TYE_pval<0.01)

pbs_gr<-makeGRangesFromDataFrame(pbs,
                                 keep.extra.columns = TRUE,
                                 seqinfo = NULL,
                                 seqnames.field=c("chromosome_name"),
                                 start.field = "start", 
                                 end.field = "end",
                                 ignore.strand = TRUE,
                                 na.rm=TRUE)
TYE_gr<-makeGRangesFromDataFrame(outliers_TYE,
                                 keep.extra.columns = TRUE,
                                 seqinfo = NULL,
                                 seqnames.field=c("chr"),
                                 start.field = "pos", 
                                 end.field = "pos",
                                 ignore.strand = TRUE,
                                 na.rm=TRUE)

TYE_sub<-subsetByOverlaps(pbs_gr, TYE_gr)
TYE_sub2<-subsetByOverlaps(TYE_gr, pbs_gr)
TYE_overlap<-as.data.frame(TYE_sub)
TYE_overlap2<-as.data.frame(TYE_sub2)
write_tsv(TYE_overlap, "TYE_pbs_sw_overlap_genes.txt")
write_tsv(TYE_overlap2, "TYE_sw_pbs_overlap_genes.txt")
TYE_overlap_pbs<-as.data.frame(cbind(TYE_overlap$seqnames, TYE_overlap$midPos, TYE_overlap$TYE_pval, TYE_overlap$ensembl_gene_id, TYE_overlap$description, TYE_overlap$external_gene_name))
names(TYE_overlap_pbs)=c("chr", "midpos", "emppval", "ensembl_geneID", "description", "gene_name")
TYE_overlap_pbs<-subset(TYE_overlap_pbs, emppval<0.01)

# Prepare the dataset
TYE.df<-na.omit(as.data.frame(TYE))
don_TYE <- TYE.df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(TYE.df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( BPcum=pos+tot) #%>%

# Add highlight_TYE and annotation information
#  mutate( is_highlight_TYE=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
#  mutate( is_annotate=ifelse(-log10(pc3_pval)>3.5, "yes", "no")) 

# Prepare X axis
axisdf <- don_TYE %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

highlight_TYE <- inner_join(don_TYE, TYE_overlap2, by = "ID")
names(highlight_TYE)=c("pos", "CLR", "alpha", "chr", "pop", "ID", "tot", "BPcum", "seqnames", "start", "end", "width", "strand", "CLR.y", "alpha.y", "pop.y")

TYE_sweep_manhat<-ggplot(don_TYE, aes(x=BPcum, y= CLR)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=0.3) +
  scale_color_manual(values = rep(c("grey", "palegreen"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  # ylim(0,8) +
  
  # Add highlight_TYEed points
  geom_point(data=highlight_TYE, color="red", size=0.3) +
  
  # Custom the theme:
  theme_bw() +
  labs(y="CLR", x="TYE") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
#TYE_sweep_manhat
#ggsave("TYE_sweep_pbs.png", TYE_sweep_manhat,width=20,height=5,dpi=600)

#################################################################################
###################### WEA ######################################################
#################################################################################

pbs<-read_tsv("../../sfs/pbs/gene_annotated_pbs.txt")
pbs<-subset(pbs, WEA_pval<0.01)

pbs_gr<-makeGRangesFromDataFrame(pbs,
                                 keep.extra.columns = TRUE,
                                 seqinfo = NULL,
                                 seqnames.field=c("chromosome_name"),
                                 start.field = "start", 
                                 end.field = "end",
                                 ignore.strand = TRUE,
                                 na.rm=TRUE)

WEA_gr<-makeGRangesFromDataFrame(outliers_WEA,
                                 keep.extra.columns = TRUE,
                                 seqinfo = NULL,
                                 seqnames.field=c("chr"),
                                 start.field = "pos", 
                                 end.field = "pos",
                                 ignore.strand = TRUE,
                                 na.rm=TRUE)

WEA_sub<-subsetByOverlaps(pbs_gr, WEA_gr)
WEA_sub2<-subsetByOverlaps(WEA_gr, pbs_gr)
WEA_overlap<-as.data.frame(WEA_sub)
WEA_overlap2<-as.data.frame(WEA_sub2)
write_tsv(WEA_overlap, "WEA_pbs_sw_overlap_genes.txt")
write_tsv(WEA_overlap2, "WEA_sw_pbs_overlap_genes.txt")
WEA_overlap_pbs<-as.data.frame(cbind(WEA_overlap$seqnames, WEA_overlap$midPos, WEA_overlap$WEA_pval, WEA_overlap$ensembl_gene_id, WEA_overlap$description, WEA_overlap$external_gene_name))
names(WEA_overlap_pbs)=c("chr", "midpos", "emppval", "ensembl_geneID", "description", "gene_name")
WEA_overlap_pbs<-subset(WEA_overlap_pbs, emppval<0.01)

# Prepare the dataset
WEA.df<-na.omit(as.data.frame(WEA))
don_WEA <- WEA.df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(WEA.df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( BPcum=pos+tot) #%>%

# Add highlight_WEA and annotation information
#  mutate( is_highlight_WEA=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
#  mutate( is_annotate=ifelse(-log10(pc3_pval)>3.5, "yes", "no")) 

# Prepare X axis
axisdf <- don_WEA %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

highlight_WEA <- inner_join(don_WEA, WEA_overlap2, by = "ID")
names(highlight_WEA)=c("pos", "CLR", "alpha", "chr", "pop", "ID", "tot", "BPcum", "seqnames", "start", "end", "width", "strand", "CLR.y", "alpha.y", "pop.y")

WEA_sweep_manhat<-ggplot(don_WEA, aes(x=BPcum, y= CLR)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=0.3) +
  scale_color_manual(values = rep(c("grey", "palegreen"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  # ylim(0,8) +
  
  # Add highlight_WEAed points
  geom_point(data=highlight_WEA, color="red", size=0.3) +
  
  # Custom the theme:
  theme_bw() +
  labs(y="CLR", x="WEA") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
#WEA_sweep_manhat
#ggsave("WEA_sweep_pbs.png", WEA_sweep_manhat,width=20,height=5,dpi=600)

#################################################################################
###################### YHE ######################################################
#################################################################################

pbs<-read_tsv("../../sfs/pbs/gene_annotated_pbs.txt")
pbs<-subset(pbs, YHE_pval<0.01)

pbs_gr<-makeGRangesFromDataFrame(pbs,
                                 keep.extra.columns = TRUE,
                                 seqinfo = NULL,
                                 seqnames.field=c("chromosome_name"),
                                 start.field = "start", 
                                 end.field = "end",
                                 ignore.strand = TRUE,
                                 na.rm=TRUE)

YHE_gr<-makeGRangesFromDataFrame(outliers_YHE,
                                 keep.extra.columns = TRUE,
                                 seqinfo = NULL,
                                 seqnames.field=c("chr"),
                                 start.field = "pos", 
                                 end.field = "pos",
                                 ignore.strand = TRUE,
                                 na.rm=TRUE)

YHE_sub<-subsetByOverlaps(pbs_gr, YHE_gr)
YHE_sub2<-subsetByOverlaps(YHE_gr, pbs_gr)
YHE_overlap<-as.data.frame(YHE_sub)
YHE_overlap2<-as.data.frame(YHE_sub2)
write_tsv(YHE_overlap, "YHE_pbs_sw_overlap_genes.txt")
write_tsv(YHE_overlap2, "YHE_sw_pbs_overlap_genes.txt")
YHE_overlap_pbs<-as.data.frame(cbind(YHE_overlap$seqnames, YHE_overlap$midPos, YHE_overlap$YHE_pval, YHE_overlap$ensembl_gene_id, YHE_overlap$description, YHE_overlap$external_gene_name))
names(YHE_overlap_pbs)=c("chr", "midpos", "emppval", "ensembl_geneID", "description", "gene_name")
YHE_overlap_pbs<-subset(YHE_overlap_pbs, emppval<0.01)

# Prepare the dataset
YHE.df<-na.omit(as.data.frame(YHE))
don_YHE <- YHE.df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(YHE.df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( BPcum=pos+tot) #%>%

# Add highlight_WEA and annotation information
#  mutate( is_highlight_WEA=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
#  mutate( is_annotate=ifelse(-log10(pc3_pval)>3.5, "yes", "no")) 

# Prepare X axis
axisdf <- don_YHE %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

highlight_YHE <- inner_join(don_YHE, YHE_overlap2, by = "ID")
names(highlight_YHE)=c("pos", "CLR", "alpha", "chr", "pop", "ID", "tot", "BPcum", "seqnames", "start", "end", "width", "strand", "CLR.y", "alpha.y", "pop.y")

YHE_sweep_manhat<-ggplot(don_YHE, aes(x=BPcum, y= CLR)) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=0.3) +
  scale_color_manual(values = rep(c("grey", "palegreen"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  # ylim(0,8) +
  
  # Add highlighted points
  geom_point(data=highlight_YHE, color="red", size=0.3) +
  
  # Custom the theme:
  theme_bw() +
  labs(y="CLR", x="YHE") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
#YHE_sweep_manhat
ggsave("YHE_sweep_pbs.png", YHE_sweep_manhat,width=20,height=5,dpi=600)

library(cowplot)

com_pbs_sweep<-plot_grid(AGA_sweep_manhat, AVB_sweep_manhat, CWB_sweep_manhat, HAY_sweep_manhat, TYD_sweep_manhat, TYE_sweep_manhat, WEA_sweep_manhat, YHE_sweep_manhat,
          ncol=1)

ggsave("com_pbs_sweep.png", com_pbs_sweep, width = 14, height =24)

