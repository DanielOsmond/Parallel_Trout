# new general run
rm(list=ls()) #clears all variables
objects() # clear all objects
graphics.off() #close all figures

#BiocManager::install("GSEABase")

setwd("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/GO")

# Load libs
lib<-c("tidyr", "dplyr","ggplot2", "data.table", "cowplot", "clusterProfiler", "enrichplot", "readr")
lapply(lib,library,character.only=T)



#GO and GSEA analysis

# Loading TERM2GENE and TERM2NAME for ALL
allg <- read.delim("term_inputs/term2gene_all.txt", header = FALSE)
alln <- read.delim("term_inputs/term2name_all.txt", header = FALSE)

# Loading TERM2GENE and TERM2NAME for BP
bpg <- read.delim("term_inputs/term2gene_BP.txt", header = FALSE)
bpn <- read.delim("term_inputs/term2name_BP.txt", header = FALSE)

# Loading TERM2GENE and TERM2NAME for MF
mfg <- read.delim("term_inputs/term2gene_MF.txt", header = FALSE)
mfn <- read.delim("term_inputs/term2name_MF.txt", header = FALSE)

# Loading TERM2GENE and TERM2NAME for CC
ccg <- read.delim("term_inputs/term2gene_CC.txt", header = FALSE)
ccn <- read.delim("term_inputs/term2name_CC.txt", header = FALSE)

### read in results for GO analysis: e.g. 
picmin02 <- read.table("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/PicMin/fdr02_genes.tsv",
                       h=T, sep = "\t", na.strings = c("", "NA"))

## perform GO enrichment: ## here the gene names are in column 6
picmin02_genes <- picmin02[,6]

picmin02_GO <- enricher(
  picmin02_genes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = ccg,
  TERM2NAME = ccn
)

picmin02_GO_cc <- as.data.frame(picmin02_GO)

picmin02_GO <- enricher(
  picmin02_genes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = bpg,
  TERM2NAME = bpn
)

picmin02_GO_bp <- as.data.frame(picmin02_GO)

picmin02_GO <- enricher(
  picmin02_genes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = mfg,
  TERM2NAME = mfn
)

picmin02_GO_mf <- as.data.frame(picmin02_GO)

write_tsv(picmin02_GO_bp, "picmin02_go_bp_fdr005.tsv")
write_tsv(picmin02_GO_mf, "picmin02_go_mf_fdr005.tsv")
write_tsv(picmin02_GO_cc, "picmin02_go_cc_fdr005.tsv")

################################################
########## Sweepfinder GO analysis #############
################################################

all_sw_pbs<-read.table("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/sweepfinder/output/merged_pbs_sw_all_genes.tsv", header=T, sep ="\t")

## perform GO enrichment:
all_sw_pbs_genes <- all_sw_pbs[,10]

all_sw_pbs_GO <- enricher(
  all_sw_pbs_genes,
  pvalueCutoff = 0.5,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = bpg,
  TERM2NAME = bpn
)

all_sw_pbs_GO_df<-as.data.frame(all_sw_pbs_GO)

### okay so not much recovered when we look at all pops - I guess that's a function of differential selection weakening power?

AGA_sw_pbs <- filter(all_sw_pbs, pop == "AGA")
AVB_sw_pbs <- filter(all_sw_pbs, pop == "AVB")
CWB_sw_pbs <- filter(all_sw_pbs, pop == "CWB")
HAY_sw_pbs <- filter(all_sw_pbs, pop == "HAY")
TYD_sw_pbs <- filter(all_sw_pbs, pop == "TYD")
TYE_sw_pbs <- filter(all_sw_pbs, pop == "TYE")
WEA_sw_pbs <- filter(all_sw_pbs, pop == "WEA")
YHE_sw_pbs <- filter(all_sw_pbs, pop == "YHE")

GO_sw_pbs <- function(df, cutoff, bpg, bpn) {
  i <- df[, 10]
  b <- enricher(
    i,
    pvalueCutoff = cutoff,
    pAdjustMethod = "fdr",
    universe = bpg,
    minGSSize = 5,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    TERM2GENE = bpg,
    TERM2NAME = bpn
  )
  
  b<-as.data.frame(b)
  return(b)
}

AGA_sw_pbs_GOb<-GO_sw_pbs(AGA_sw_pbs, 0.1, bpg, bpn)
AVB_sw_pbs_GOb<-GO_sw_pbs(AVB_sw_pbs, 0.1, bpg, bpn)
CWB_sw_pbs_GOb<-GO_sw_pbs(CWB_sw_pbs, 0.1, bpg, bpn)
HAY_sw_pbs_GOb<-GO_sw_pbs(HAY_sw_pbs, 0.1, bpg, bpn)
TYD_sw_pbs_GOb<-GO_sw_pbs(TYD_sw_pbs, 0.1, bpg, bpn)
TYE_sw_pbs_GOb<-GO_sw_pbs(TYE_sw_pbs, 0.1, bpg, bpn)
WEA_sw_pbs_GOb<-GO_sw_pbs(WEA_sw_pbs, 0.1, bpg, bpn)
YHE_sw_pbs_GOb<-GO_sw_pbs(YHE_sw_pbs, 0.1, bpg, bpn)

pop_names <- c("AGA", "AVB", "CWB", "HAY", "TYD", "TYE", "WEA", "YHE")
dataframes <- lapply(pop_names, function(pop) get(paste0(pop, "_sw_pbs_GOb")))

# Add a 'pop' column to each dataframe
dataframes_with_pop <- Map(cbind, dataframes, pop = pop_names)

# Combine dataframes into a single dataframe
merged_sw_pbs_GOb <- bind_rows(dataframes_with_pop)
merged_sw_pbs_GOb = data.table(merged_sw_pbs_GO)[order(p.adjust),]

AGA_sw_pbs_GOm<-GO_sw_pbs(AGA_sw_pbs, 0.1, mfg, mfn)
AVB_sw_pbs_GOm<-GO_sw_pbs(AVB_sw_pbs, 0.1, mfg, mfn)
CWB_sw_pbs_GOm<-GO_sw_pbs(CWB_sw_pbs, 0.1, mfg, mfn)
HAY_sw_pbs_GOm<-GO_sw_pbs(HAY_sw_pbs, 0.1, mfg, mfn)
TYD_sw_pbs_GOm<-GO_sw_pbs(TYD_sw_pbs, 0.1, mfg, mfn)
TYE_sw_pbs_GOm<-GO_sw_pbs(TYE_sw_pbs, 0.1, mfg, mfn)
WEA_sw_pbs_GOm<-GO_sw_pbs(WEA_sw_pbs, 0.1, mfg, mfn)
YHE_sw_pbs_GOm<-GO_sw_pbs(YHE_sw_pbs, 0.1, mfg, mfn)

pop_names <- c("AGA", "AVB", "CWB", "HAY", "TYD", "TYE", "WEA", "YHE")
dataframes <- lapply(pop_names, function(pop) get(paste0(pop, "_sw_pbs_GOm")))

# Add a 'pop' column to each dataframe
dataframes_with_pop <- Map(cbind, dataframes, pop = pop_names)

# Combine dataframes into a single dataframe
merged_sw_pbs_GOm <- bind_rows(dataframes_with_pop)
merged_sw_pbs_GOm = data.table(merged_sw_pbs_GOm)[order(p.adjust),]


AGA_sw_pbs_GOc<-GO_sw_pbs(AGA_sw_pbs, 0.1, ccg, ccn)
AVB_sw_pbs_GOc<-GO_sw_pbs(AVB_sw_pbs, 0.1, ccg, ccn)
CWB_sw_pbs_GOc<-GO_sw_pbs(CWB_sw_pbs, 0.1, ccg, ccn)
HAY_sw_pbs_GOc<-GO_sw_pbs(HAY_sw_pbs, 0.1, ccg, ccn)
TYD_sw_pbs_GOc<-GO_sw_pbs(TYD_sw_pbs, 0.1, ccg, ccn)
TYE_sw_pbs_GOc<-GO_sw_pbs(TYE_sw_pbs, 0.1, ccg, ccn)
WEA_sw_pbs_GOc<-GO_sw_pbs(WEA_sw_pbs, 0.1, ccg, ccn)
YHE_sw_pbs_GOc<-GO_sw_pbs(YHE_sw_pbs, 0.1, ccg, ccn)

pop_names <- c("AGA", "AVB", "CWB", "HAY", "TYE", "WEA", "YHE")
dataframes <- lapply(pop_names, function(pop) get(paste0(pop, "_sw_pbs_GOc")))

# Add a 'pop' column to each dataframe
dataframes_with_pop <- Map(cbind, dataframes, pop = pop_names)

# Combine dataframes into a single dataframe
merged_sw_pbs_GOc <- bind_rows(dataframes_with_pop)
merged_sw_pbs_GOc = data.table(merged_sw_pbs_GOc)[order(p.adjust),]

write_tsv(merged_sw_pbs_GOb, "merged_sw_pbs_GO_bp.tsv")
write_tsv(merged_sw_pbs_GOm, "merged_sw_pbs_GO_mf.tsv")
write_tsv(merged_sw_pbs_GOc, "merged_sw_pbs_GO_cc.tsv")


## AF_pbs intersect ## 

Wales_AF_pbs <- read.table("../cymru.AF_fully_pbs_overlap.tsv", header = T, sep = "\t", na.strings = c("", "NA"))

Wales_AF_pbs_genes <- Wales_AF_pbs[,17]

Wales_AF_pbs_bp <- enricher(
  Wales_AF_pbs_genes,
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = mfg,
  TERM2NAME = mfn
)

Wales_AF_pbs_bp <- as.data.frame(Wales_AF_pbs_bp)

###################################################################################################################
### AF-vapeR ######################################################################################################
###################################################################################################################

### UK ###

UK_AF<- read.table("../mafs/UK_eig1_99_fully_genes.tsv", header = T, sep = "\t", na.strings = c("", "NA"))

## perform GO enrichment: ## here the gene names are in column 1
UK_AF_genes <- UK_AF[,6]

UK_AF_bp <- enricher(
  UK_AF_genes,
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = bpg,
  TERM2NAME = bpn
)

UK_AF_bp <- as.data.frame(UK_AF_bp)

UK_AF_mf <- enricher(
  UK_AF_genes,
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = mfg,
  TERM2NAME = mfn
)

UK_AF_mf <- as.data.frame(UK_AF_mf)

UK_AF_cc <- enricher(
  UK_AF_genes,
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = ccg,
  TERM2NAME = ccn
)

UK_AF_cc <- as.data.frame(UK_AF_cc)

write_tsv(UK_AF_bp, "UK_AF_bp.tsv")
write_tsv(UK_AF_mf, "UK_AF_mf.tsv")
write_tsv(UK_AF_cc, "UK_AF_cc.tsv")

### Wales ###

Wales_AF<- read.table("../mafs/cymru_eig1_99_fully_genes.tsv", header = T, sep = "\t", na.strings = c("", "NA"))

## perform GO enrichment: ## here the gene names are in column 1
Wales_AF_genes <- Wales_AF[,6]

Wales_AF_bp <- enricher(
  Wales_AF_genes,
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = bpg,
  TERM2NAME = bpn
)

Wales_AF_bp <- as.data.frame(Wales_AF_bp)

Wales_AF_mf <- enricher(
  Wales_AF_genes,
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = mfg,
  TERM2NAME = mfn
)

Wales_AF_mf <- as.data.frame(Wales_AF_mf)

Wales_AF_cc <- enricher(
  Wales_AF_genes,
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = ccg,
  TERM2NAME = ccn
)

Wales_AF_cc <- as.data.frame(Wales_AF_cc)

write_tsv(Wales_AF_bp, "Wales_AF_bp.tsv")
write_tsv(Wales_AF_mf, "Wales_AF_mf.tsv")
write_tsv(Wales_AF_cc, "Wales_AF_cc.tsv")

### not much really going on here with regards to metals, mostly oxygen transport / immune function.

### NE England ###

NE_AF<- read.table("../mafs/NE_eig1_99_fully_genes.tsv", header = T, sep = "\t", na.strings = c("", "NA"))

## perform GO enrichment: ## here the gene names are in column 1
NE_AF_genes <- NE_AF[,6]

NE_AF_bp <- enricher(
  NE_AF_genes,
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = bpg,
  TERM2NAME = bpn
)

NE_AF_bp <- as.data.frame(NE_AF_bp)

NE_AF_mf <- enricher(
  NE_AF_genes,
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = mfg,
  TERM2NAME = mfn
)

NE_AF_mf <- as.data.frame(NE_AF_mf)

NE_AF_cc <- enricher(
  NE_AF_genes,
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = ccg,
  TERM2NAME = ccn
)

NE_AF_cc <- as.data.frame(NE_AF_cc)

write_tsv(NE_AF_bp, "NE_AF_bp.tsv")
write_tsv(NE_AF_mf, "NE_AF_mf.tsv")
write_tsv(NE_AF_cc, "NE_AF_cc.tsv")

#########################################
##### Water Chem ########################
#########################################

### Cu ###

Cu_AF<- read.table("../mafs/Cu_eig1_99_fully_genes.tsv", header = T, sep = "\t", na.strings = c("", "NA"))

## perform GO enrichment: ## here the gene names are in column 1
Cu_AF_genes <- Cu_AF[,6]

Cu_AF_bp <- enricher(
  Cu_AF_genes,
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = bpg,
  TERM2NAME = bpn
)

Cu_AF_bp <- as.data.frame(Cu_AF_bp)

Cu_AF_mf <- enricher(
  Cu_AF_genes,
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = mfg,
  TERM2NAME = mfn
)

Cu_AF_mf <- as.data.frame(Cu_AF_mf)

Cu_AF_cc <- enricher(
  Cu_AF_genes,
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = ccg,
  TERM2NAME = ccn
)

Cu_AF_cc <- as.data.frame(Cu_AF_cc)

write_tsv(Cu_AF_bp, "Cu_AF_bp.tsv")
write_tsv(Cu_AF_mf, "Cu_AF_mf.tsv")
write_tsv(Cu_AF_cc, "Cu_AF_cc.tsv")


### Pb ###

Pb_AF<- read.table("../mafs/Pb_eig1_99_fully_genes.tsv", header = T, sep = "\t", na.strings = c("", "NA"))

## perform GO enrichment: ## here the gene names are in column 1
Pb_AF_genes <- Pb_AF[,6]

Pb_AF_bp <- enricher(
  Pb_AF_genes,
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = bpg,
  TERM2NAME = bpn
)

Pb_AF_bp <- as.data.frame(Pb_AF_bp)

Pb_AF_mf <- enricher(
  Pb_AF_genes,
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = mfg,
  TERM2NAME = mfn
)

Pb_AF_mf <- as.data.frame(Pb_AF_mf)

Pb_AF_cc <- enricher(
  Pb_AF_genes,
  pvalueCutoff = 0.1,
  pAdjustMethod = "fdr",
  universe=bpg,
  minGSSize = 5,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  TERM2GENE = ccg,
  TERM2NAME = ccn
)

Pb_AF_cc <- as.data.frame(Pb_AF_cc)

write_tsv(Pb_AF_bp, "Pb_AF_bp.tsv")
write_tsv(Pb_AF_mf, "Pb_AF_mf.tsv")
write_tsv(Pb_AF_cc, "Pb_AF_cc.tsv")
