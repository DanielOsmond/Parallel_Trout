### PicMin script provided by Jim Whiting Feb 2023 ### 

# An example PicMin analysis...
library(poolr)
library(ggplot2)
library(data.table)
library(readr)
setwd("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/PicMin")
source("./PicMin.R")
library("qvalue")
library(ape)
library(GenomicRanges)
library(Rfast)


# read in the Fst values per pop and calculate empPvals 
# skip this for 100k window set where run on isca

# anglesey 
r<-read.delim("50k_out/AGA_AGB_win50k_subset.fst.txt",sep="\t",head=F)
names(r)[c(1:5)] <- c("comp","region","midPos","Nsites","uwFst")
r$pval1 <- qvalue::empPvals(stat = r$uwFst,stat0 = r$uwFst)
hist (r$pval) #should be a uniform distribution of p - just ranking values of Fst rather than wanting stat signif
head(r)
# write new file  #
write_tsv(r, "AGA_AGB_50k_fst.pvals.txt")

# Ireland1 
r<-read.delim("50k_out/AVB_AVD_win50k_subset.fst.txt",sep="\t",head=F)
names(r)[c(1:5)] <- c("comp","region","midPos","Nsites","uwFst")
r$pval1 <- qvalue::empPvals(stat = r$uwFst,stat0 = r$uwFst)
hist (r$pval) #should be a uniform distribution of p - just ranking values of Fst rather than wanting stat signif
head(r)
# write new file  #
write_tsv(r, "AVB_AVD_50k_fst.pvals.txt")

#Ireland2
r<-read.delim("50k_out/AVB_OWE_win50k_subset.fst.txt",sep="\t",head=F)
names(r)[c(1:5)] <- c("comp","region","midPos","Nsites","uwFst")
r$pval1 <- qvalue::empPvals(stat = r$uwFst,stat0 = r$uwFst)
hist (r$pval) #should be a uniform distribution of p - just ranking values of Fst rather than wanting stat signif
head(r)
# write new file  #
write_tsv(r, "AVB_OWE_50k_fst.pvals.txt")

#NWales1
r<-read.delim("50k_out/CWA_CWC_win50k_subset.fst.txt",sep="\t",head=F)
names(r)[c(1:5)] <- c("comp","region","midPos","Nsites","uwFst")
r$pval1 <- qvalue::empPvals(stat = r$uwFst,stat0 = r$uwFst)
hist (r$pval) #should be a uniform distribution of p - just ranking values of Fst rather than wanting stat signif
head(r)
# write new file  #
write_tsv(r, "CWA_CWC_50k_fst.pvals.txt")

#NWales2
r<-read.delim("50k_out/CWB_CWC_win50k_subset.fst.txt",sep="\t",head=F)
names(r)[c(1:5)] <- c("comp","region","midPos","Nsites","uwFst")
r$pval1 <- qvalue::empPvals(stat = r$uwFst,stat0 = r$uwFst)
hist (r$pval) #should be a uniform distribution of p - just ranking values of Fst rather than wanting stat signif
head(r)
# write new file  #
write_tsv(r, "CWB_CWC_50k_fst.pvals.txt")

#Cornwall
r<-read.delim("50k_out/HAY1_CAM_win50k_subset.fst.txt",sep="\t",head=F)
names(r)[c(1:5)] <- c("comp","region","midPos","Nsites","uwFst")
r$pval1 <- qvalue::empPvals(stat = r$uwFst,stat0 = r$uwFst)
hist (r$pval) #should be a uniform distribution of p - just ranking values of Fst rather than wanting stat signif
head(r)
# write new file  #
write_tsv(r, "HAY1_CAM_50k_fst.pvals.txt")

#Tyne1
r<-read.delim("50k_out/TYA_TYD_win50k_subset.fst.txt",sep="\t",head=F)
names(r)[c(1:5)] <- c("comp","region","midPos","Nsites","uwFst")
r$pval1 <- qvalue::empPvals(stat = r$uwFst,stat0 = r$uwFst)
hist (r$pval) #should be a uniform distribution of p - just ranking values of Fst rather than wanting stat signif
head(r)
# write new file  #
write_tsv(r, "TYA_TYD_50k_fst.pvals.txt")

#Tyne2
r<-read.delim("50k_out/TYA_TYE_win50k_subset.fst.txt",sep="\t",head=F)
names(r)[c(1:5)] <- c("comp","region","midPos","Nsites","uwFst")
r$pval1 <- qvalue::empPvals(stat = r$uwFst,stat0 = r$uwFst)
hist (r$pval) #should be a uniform distribution of p - just ranking values of Fst rather than wanting stat signif
head(r)
# write new file  #
write_tsv(r, "TYA_TYE_50k_fst.pvals.txt")

#Wear
r<-read.delim("50k_out/WEA_WEF_win50k_subset.fst.txt",sep="\t",head=F)
names(r)[c(1:5)] <- c("comp","region","midPos","Nsites","uwFst")
r$pval1 <- qvalue::empPvals(stat = r$uwFst,stat0 = r$uwFst)
hist (r$pval) #should be a uniform distribution of p - just ranking values of Fst rather than wanting stat signif
head(r)
# write new file  #
write_tsv(r, "WEA_WEF_50k_fst.pvals.txt")

#Mid_Wales
r<-read.delim("50k_out/YHE_NEV2_win50k_subset.fst.txt",sep="\t",head=F)
names(r)[c(1:5)] <- c("comp","region","midPos","Nsites","uwFst")
r$pval1 <- qvalue::empPvals(stat = r$uwFst,stat0 = r$uwFst)
hist (r$pval) #should be a uniform distribution of p - just ranking values of Fst rather than wanting stat signif
head(r)
# write new file  #
write_tsv(r, "YHE_NEV2_50k_fst.pvals.txt")

# Read in the Fst Pvals per population 
Anglesey<- read_tsv("AGA_AGB_50k_fst.pvals.txt")
Ireland1<- read_tsv("AVB_AVD_50k_fst.pvals.txt")
Ireland2<- read_tsv("AVB_OWE_50k_fst.pvals.txt")
NWales1<-  read_tsv("CWA_CWC_50k_fst.pvals.txt")
NWales2<-  read_tsv("CWB_CWC_50k_fst.pvals.txt")
Kernow<-   read_tsv("HAY1_CAM_50k_fst.pvals.txt")
Tyne1<-    read_tsv("TYA_TYD_50k_fst.pvals.txt")
Tyne2<-    read_tsv("TYA_TYE_50k_fst.pvals.txt")
Wear<-     read_tsv("WEA_WEF_50k_fst.pvals.txt")
Mid_Wales<-read_tsv("YHE_NEV2_50k_fst.pvals.txt")

## gff filtering 
## filtering Fst windows that have genes intersecting them to reduce number of comparisons and thereby FDR correction

#early stuff - ignore #
#trout_gff<-read.gff("Salmo_trutta.fSalTru1.1.103.gff3", na.strings = c(".", "?"), GFF3 = TRUE)
#trout_gff$chr<-as.numeric(trout_gff$seqid)
#chr_gff<-data.frame(subset(trout_gff, chr>=0| chr<=40, type=-c("region")))

## naughty here and manually entered header in notepad 
## created the bed using gff2bed from BEDOPs on isca
gff_bed<-read_tsv("trout_gff.bed")
#names(gff_bed) =c("v1","v2","v3","v4","v5","v6","v7","cat","v9","v10")

#have a wee look - seems like stuff all read in as character rather than factor
summary(gff_bed)
#gff_bed<-as.data.frame(gff_bed)
gff_bed$type<-as.factor(gff_bed$type)
gff_bed$chr<-as.factor(gff_bed$chr)
gff_bed$v7<-as.factor(gff_bed$v7)

#get rid of things we're not interested in

gff_bed_genes <- data.frame(subset(gff_bed, type!="region"))
gff_bed_genes <-data.frame(subset(gff_bed_genes, type!="biological_region"))
gff_bed_genes <-data.frame(subset(gff_bed, type=="gene"))

#check removal
summary(gff_bed_genes$v7)

#look at the size of genes remaining
gff_bed_genes$len<-gff_bed_genes$end - gff_bed_genes$start
summary(gff_bed_genes$len)

gene_ranges = subset(gff_bed_genes, select= c(1,2,3))
write_tsv(gene_ranges, "gene_ranges.tsv")


#just use Anglesey as proxy for all windowed Fst in that window range. 
Anglesey<-read.delim("50k_out/AGA_AGB_win50k_subset.fst.txt",sep="\t",head=F)
names(Anglesey)[c(1:5)] <- c("comp","region","midPos","Nsites","uwFst")

Anglesey$start <- Anglesey$midPos - 25000
Anglesey$end <- Anglesey$midPos + 25000

gr1<-makeGRangesFromDataFrame(gene_ranges,
                         keep.extra.columns = FALSE,
                         seqinfo = NULL,
                         seqnames.field=c("chr"),
                         start.field = "start", 
                         end.field = "end",
                         ignore.strand = TRUE,
                         na.rm=TRUE)

gr2<-makeGRangesFromDataFrame(Anglesey,
                              keep.extra.columns = FALSE,
                              seqinfo = NULL,
                              seqnames.field=c("region"),
                              start.field = "start", 
                              end.field = "end",
                              ignore.strand = TRUE,
                              na.rm=TRUE)

gr_sub<-subsetByOverlaps(gr2,gr1)
gr_dat_all_50k<-as.data.frame(gr_sub)
write_tsv(gr_dat_all_50k, "gr_dat_all_50k.tsv")

summary(gr_dat$seqnames)

Anglesey$window<-paste0(Anglesey$region,":",Anglesey$midPos)

### read back in windows containing ensembl annotated genes ### 

windows<-read.table("gr_dat_all_50k.tsv", header=T, sep = "\t")
windows$midPos <- (windows$start + 25000)

#combine these to a dataframe
# use same vector sets as AF-VapeR
#uk wide 
## read in as uwFst matrix and then recalc empPvals post filtering ## 
uwfst_res<-cbind.data.frame(Anglesey[2], Anglesey[3], 
                          Anglesey[5], Ireland1[5], 
                          Ireland2[5],
                          NWales1[5], NWales2[5], 
                          Kernow[5],
                          Tyne1[5], Tyne2[5], 
                          Wear[5], Mid_Wales[5])

#uk labels
names(uwfst_res) = c("chr", "midPos", 
                    "Anglesey", "Ireland1", 
                    "Ireland2", "NWales1",
                   "NWales2",
                   "Kernow", "Tyne1", "Tyne2", 
                   "Wear", 
                   "Mid_Wales")
uwfst_res$ID <- paste0(uwfst_res$chr, "_", uwfst_res$midPos)
windows$ID <- paste0(windows$seqnames, "_", windows$midPos)

fst_res<-setDT(uwfst_res)[ID %chin% windows$ID]
uwfst_genes<-as.data.frame(fst_res)

write_tsv(uwfst_genes, "uwFst_matrix_50k_genes.tsv")

uwfst_genes$Anglesey_p <- qvalue::empPvals(stat = uwfst_genes$Anglesey,stat0 = uwfst_genes$Anglesey)
uwfst_genes$Ireland1_p <- qvalue::empPvals(stat = uwfst_genes$Ireland1,stat0 = uwfst_genes$Ireland1)
uwfst_genes$Ireland2_p <- qvalue::empPvals(stat = uwfst_genes$Ireland2,stat0 = uwfst_genes$Ireland2)
uwfst_genes$Nwales1_p <- qvalue::empPvals(stat = uwfst_genes$NWales1,stat0 = uwfst_genes$NWales1)
uwfst_genes$NWales2_p <- qvalue::empPvals(stat = uwfst_genes$NWales2,stat0 = uwfst_genes$NWales2)
uwfst_genes$Kernow_p <- qvalue::empPvals(stat = uwfst_genes$Kernow,stat0 = uwfst_genes$Kernow)
uwfst_genes$Tyne1_p <- qvalue::empPvals(stat = uwfst_genes$Tyne1,stat0 = uwfst_genes$Tyne1)
uwfst_genes$Tyne2_p <- qvalue::empPvals(stat = uwfst_genes$Tyne2,stat0 = uwfst_genes$Tyne2)
uwfst_genes$Wear_p <- qvalue::empPvals(stat = uwfst_genes$Wear,stat0 = uwfst_genes$Wear)
uwfst_genes$MidWales_p <- qvalue::empPvals(stat = uwfst_genes$Mid_Wales,stat0 = uwfst_genes$Mid_Wales)

empP_pops<- uwfst_genes[,13:23]

write_tsv(empP_pops, "empPvals_matrix_50k_genes.tsv")
### now read in the columns that we actually want! ### 

fst_mat<-read_tsv("empPvals_matrix_50k_genes.tsv")

#for UK_wide 7 pairs#
fst_res<-as.data.frame(subset(fst_mat, select =-c(3,6,9)))

#all_pairs
fst_res<-as.data.frame(fst_mat)


# Consider that this matrix contains 7 Fst analyses for 10,000 regions of the genome
# We have converted each of the columns into a vector of empirical pvals
# So each column is a perfectly uniform distribution
#fst_res = matrix(ncol = 7,nrow = 10000)
#for(col in 1:ncol(fst_res)){
#  fst_res[,col] = runif(nrow(fst_res))
#}
# For e.g. the empirical pvals from 1 population pair MUST look like this
# basically a flat distribution
hist(fst_res[,2])

### version 2 ###

# Running PicMin ----------------------------------------------------------
# First, we need to define the correlation structure of ordered pvals
# This is because when drawing 7 random pvals and ordering them, there
# will be a correlation structure among the pvals that needs accounting for

### The reason we needed the correlation matrix was for estimating the tippett empirical null distribution. But actually we can just do this once and save loads of time and effort.
### New function for doing this also does the correlation matrix step for us
### This function now also includes a seed for repeatability across analyses
seedN = 1966
tippett_nulls = generateTippettNull(N = 7,
                                    geneN = 40000,
                                    seed = seedN,
                                    numReps = 1000000)

### Running PicMin is now much quicker and straightforward, we just need to give a vector of pvals and a vector of 'groups' in your case window ids
fst_input = melt(fst_res)
colnames(fst_input) = c('window','replicate','epvalue')
picmin_res = PicMin_bulk(pvals = fst_input$epvalue,
                         groups = fst_input$window,
                         nulls = tippett_nulls)

# This output then gives a pvalue for the picmin test per window, and the 'config.est' gives 
# an approximation of how many of the 7 tests are 'contributing' towards the repeatability here
# Specifically, the config.est is the pvalue that is the 'least expected' relative to the expectation
# of the Nth pvalue.

# The first inspection of this is which are strongest
p_ordered<-picmin_res[order(picmin_p),]
picmin_res[2890:2910,]
# Bare in mind that we are doing LOTS of tests, so will find some very low pvalues
# This means that we should do an FDR correction really
picmin_res$fdr = p.adjust(picmin_res$p,method = 'fdr')
hist(picmin_res$fdr, xlim=c(0,0.2), ylim =c(0,100)) # Nothing significant here, good because the data is all random!

write_tsv(p_ordered, "uk_50k_picmin_pordered.tsv")
# Be aware with the fdr correction, that the test will lose power the more windows in the genome
# that you look at. Because of this, the most desirable thing to do is minimise the number of tests
# Maybe consider testing larger windows, 50kb instead of 10kb, and maybe only look at windows overlapping
# genes if that's what you're really interested in!!

# Can also visualise these results in a manhattan plot of course which might be interesting
picmin_res$index = 1:nrow(picmin_res)
ggplot(picmin_res,aes(index,-log10(picmin_p))) +
  geom_point()

##### The below is still valid btw for testing the Nth pval, and doesn't require a tippett correction so don't worry
# Alternative PicMin ------------------------------------------------------
# An alternative approach is to only test the Nth pval, and not test pvals 2-N
# This is more powerful because we don't have to multiple correct for the N-1 tests
# However it tests a specific all or nothing hypothesis, like af-vaper, i.e. ALL populations
# are differentiating the same locus. So it's a trade-off of hypothesis flexibility and
# statistical power. It's also a lot quicker to run because we don't have to multiple correct 
# And do all those internal perms!
alternative_picmin_res <- rbindlist(lapply(1:nrow(fst_res[,2:8]),function(i){
  print(i)
  set.seed(1000)
  
  # We don't actually need the pval corr mat or numreps here because they are ignored
  # We fix the config estimate at the maximum value to specifically test the 7th pval
  PicMin(fst_res[i,], pval_corr_mat, numReps = 1000,fixConfig = ncol(fst_res))
}))

# Can treat these result as before, but be aware the hypothesis being tested is different!!
alternative_picmin_res[order(p),]

# Be aware with the fdr correction, that the test will lose power the more windows in the genome
# that you look at. Because of this, the most desirable thing to do is minimise the number of tests
# Maybe consider testing larger windows, 50kb instead of 10kb, and maybe only look at windows overlapping
# genes if that's what you're really interested in!!

# Can also visualise these results in a manhattan plot of course which might be interesting
picmin_res$index = 1:nrow(picmin_res)
ggplot(picmin_res,aes(index,-log10(p))) +
  geom_point()


### look at Fst regions of interest ### 
setwd("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/PicMin/per_locus")

AGA_AGB<-read.table("AGA_AGB.alpha_beta.txt", header=F, sep="\t")
AVB_AVD<-read.table("AVB_AVD.alpha_beta.txt", header=F, sep="\t")
AVB_OWE<-read.table("AVB_OWE.alpha_beta.txt", header=F, sep="\t")
CWA_CWC<-read.table("CWA_CWC.alpha_beta.txt", header=F, sep="\t")
CWB_CWC<-read.table("CWB_CWC.alpha_beta.txt", header=F, sep="\t")
HAY_CAM<-read.table("HAY1_CAM.alpha_beta.txt", header=F, sep="\t")
TYA_TYD<-read.table("TYA_TYD.alpha_beta.txt", header=F, sep="\t")
TYA_TYE<-read.table("TYA_TYE.alpha_beta.txt", header=F, sep="\t")
WEA_WEF<-read.table("WEA_WEF.alpha_beta.txt", header=F, sep="\t")
YHE_NEV<-read.table("YHE_NEV2.alpha_beta.txt", header=F, sep="\t")

names(AGA_AGB)=c("chr", "pos", "uWFst", "wFst")
names(AVB_AVD)=c("chr", "pos", "uWFst", "wFst")
names(AVB_OWE)=c("chr", "pos", "uWFst", "wFst")
names(CWA_CWC)=c("chr", "pos", "uWFst", "wFst")
names(CWB_CWC)=c("chr", "pos", "uWFst", "wFst")
names(HAY_CAM)=c("chr", "pos", "uWFst", "wFst")
names(TYA_TYD)=c("chr", "pos", "uWFst", "wFst")
names(TYA_TYE)=c("chr", "pos", "uWFst", "wFst")
names(WEA_WEF)=c("chr", "pos", "uWFst", "wFst")
names(YHE_NEV)=c("chr", "pos", "uWFst", "wFst")

library("qvalue")
library("readr")
TYA_TYD$pval1 <- qvalue::empPvals(stat = TYA_TYD$uWFst,stat0 = TYA_TYD$uWFst)
AGA_AGB$pval1 <- qvalue::empPvals(stat = AGA_AGB$uWFst,stat0 = AGA_AGB$uWFst)
WEA_WEF$pval1 <- qvalue::empPvals(stat = WEA_WEF$uWFst,stat0 = WEA_WEF$uWFst)
hist(TYA_TYD$pval1)

AGA_AGB_81<- data.frame(subset(AGA_AGB, chr =="NC_042981.1"))
TYA_TYD_81<- data.frame(subset(TYA_TYD, chr =="NC_042981.1"))
WEA_WEF_81<- data.frame(subset(WEA_WEF, chr =="NC_042981.1"))
#  pos > 23550000 & pos <25750000

library(ggplot2)
source("/Users/do271/OneDrive - University of Exeter/PhD/Genomics/pcangsd/manhatten_plot.R")

anglesey_81_manhat<- ggplot(AGA_AGB_81, aes(x=pos, y=uWFst) +
                              geom_point(alpha = 0.75),
                            xlim =c(23550000, 25750000))
plot(AGA_AGB_81$uWFst~AGA_AGB_81$pos, xlim=c(23550000, 25750000))                            
plot(TYA_TYD_81$uWFst~TYA_TYD_81$pos, xlim=c(23550000, 25750000))
par(mfrow = c(3,1))
plot((-log10(AGA_AGB_81$pval1))~AGA_AGB_81$pos, xlim=c(23550000, 25750000), ylab= "-log10(p)")
plot((-log10(TYA_TYD_81$pval1))~TYA_TYD_81$pos, xlim=c(23550000, 25750000), ylab= "-log10(p)")
plot((-log10(WEA_WEF_81$pval1))~WEA_WEF_81$pos, xlim=c(23550000, 25750000), ylab= "-log10(p)")

#### PBS Analysis ####
#### Calculate from tutorial #### 

# read in the files 
setwd("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/sfs/pbs")

library(ggplot2)
library(qvalue)
library(readr)
library(GenomicRanges)
library(data.table)

AGA_pbs<-read.table("AGA_chr.pbs.fst.txt", sep = '\t', header=F)
AVB_pbs<-read.table("AVB_chr.pbs.fst.txt", sep = '\t', header=F)
HAY_pbs<-read.table("HAY1_chr.pbs.fst.txt", sep = '\t', header=F)
CWB_pbs<-read.table("CWB_chr.pbs.fst.txt", sep = '\t', header=F)
TYD_pbs<-read.table("TYD_chr.pbs.fst.txt", sep = '\t', header=F)
TYE_pbs<-read.table("TYE_chr.pbs.fst.txt", sep = '\t', header=F)
WEA_pbs<-read.table("WEA_chr.pbs.fst.txt", sep = '\t', header=F)
YHE_pbs<-read.table("YHE_chr.pbs.fst.txt", sep = '\t', header=F)

names(AGA_pbs) = c("region",  "chr", "midPos",  "Nsites",  "Fst01",   "Fst02",   "Fst12", "PBS0", "PBS1", "PBS2")
AGA_pbs$pval<-qvalue::empPvals(stat = AGA_pbs$PBS0,stat0 = AGA_pbs$PBS0)

names(AVB_pbs) = c("region",  "chr", "midPos",  "Nsites",  "Fst01",   "Fst02",   "Fst12", "PBS0", "PBS1", "PBS2")
AVB_pbs$pval<-qvalue::empPvals(stat = AVB_pbs$PBS0,stat0 = AVB_pbs$PBS0)

names(HAY_pbs) = c("region",  "chr", "midPos",  "Nsites",  "Fst01",   "Fst02",   "Fst12", "PBS0", "PBS1", "PBS2")
HAY_pbs$pval<-qvalue::empPvals(stat = HAY_pbs$PBS0,stat0 = HAY_pbs$PBS0)

names(CWB_pbs) = c("region",  "chr", "midPos",  "Nsites",  "Fst01",   "Fst02",   "Fst12", "PBS0", "PBS1", "PBS2")
CWB_pbs$pval<-qvalue::empPvals(stat = CWB_pbs$PBS0,stat0 = CWB_pbs$PBS0)

names(TYD_pbs) = c("region",  "chr", "midPos",  "Nsites",  "Fst01",   "Fst02",   "Fst12", "PBS0", "PBS1", "PBS2")
TYD_pbs$pval<-qvalue::empPvals(stat = TYD_pbs$PBS0,stat0 = TYD_pbs$PBS0)

names(TYE_pbs) = c("region",  "chr", "midPos",  "Nsites",  "Fst01",   "Fst02",   "Fst12", "PBS0", "PBS1", "PBS2")
TYE_pbs$pval<-qvalue::empPvals(stat = TYE_pbs$PBS0,stat0 = TYE_pbs$PBS0)

names(WEA_pbs) = c("region",  "chr", "midPos",  "Nsites",  "Fst01",   "Fst02",   "Fst12", "PBS0", "PBS1", "PBS2")
WEA_pbs$pval<-qvalue::empPvals(stat = WEA_pbs$PBS0,stat0 = WEA_pbs$PBS0)

names(YHE_pbs) = c("region",  "chr", "midPos",  "Nsites",  "Fst01",   "Fst02",   "Fst12", "PBS0", "PBS1", "PBS2")
YHE_pbs$pval<-qvalue::empPvals(stat = YHE_pbs$PBS0,stat0 = YHE_pbs$PBS0)

#plot those puppies# 

AGA_pbs_manhat<-ggplot(AGA_pbs,aes(-log10(pval),x = midPos)) +
  geom_point() +
  facet_wrap(~chr,scales = 'free_x')
ggsave("AGA_pbs_manhat.png", plot=AGA_pbs_manhat, width =30, height = 20)

AVB_pbs_manhat<-ggplot(AVB_pbs,aes(-log10(pval),x = midPos)) +
  geom_point() +
  facet_wrap(~chr,scales = 'free_x')
ggsave("AVB_pbs_manhat.png", plot=AVB_pbs_manhat, width =30, height = 20)

HAY_pbs_manhat<-ggplot(HAY_pbs,aes(-log10(pval),x = midPos)) +
  geom_point() +
  facet_wrap(~chr,scales = 'free_x')
ggsave("HAY_pbs_manhat.png", plot=HAY_pbs_manhat, width =30, height = 20)

CWB_pbs_manhat<-ggplot(CWB_pbs,aes(-log10(pval),x = midPos)) +
  geom_point() +
  facet_wrap(~chr,scales = 'free_x')
ggsave("CWB_pbs_manhat.png", plot=CWB_pbs_manhat, width =30, height = 20)

TYD_pbs_manhat<-ggplot(TYD_pbs,aes(-log10(pval),x = midPos)) +
  geom_point() +
  facet_wrap(~chr,scales = 'free_x')
ggsave("TYD_pbs_manhat.png", plot=TYD_pbs_manhat, width =30, height = 20)

TYE_pbs_manhat<-ggplot(TYE_pbs,aes(-log10(pval),x = midPos)) +
  geom_point() +
  facet_wrap(~chr,scales = 'free_x')
ggsave("TYE_pbs_manhat.png", plot=TYE_pbs_manhat, width =30, height = 20)

WEA_pbs_manhat<-ggplot(WEA_pbs,aes(-log10(pval),x = midPos)) +
  geom_point() +
  facet_wrap(~chr,scales = 'free_x')
ggsave("WEA_pbs_manhat.png", plot=WEA_pbs_manhat, width =30, height = 20)

YHE_pbs_manhat<-ggplot(YHE_pbs,aes(-log10(pval),x = midPos)) +
  geom_point() +
  facet_wrap(~chr,scales = 'free_x')
ggsave("YHE_pbs_manhat.png", plot=YHE_pbs_manhat, width =30, height = 20)


# combine these into a matrix with a pbs stat per metal pop 
pbs_mat<- as.data.frame(cbind(YHE_pbs$chr, YHE_pbs$midPos, AGA_pbs$PBS0, AVB_pbs$PBS0, CWB_pbs$PBS0, HAY_pbs$PBS0,
               TYD_pbs$PBS0, TYE_pbs$PBS0, WEA_pbs$PBS0, YHE_pbs$PBS0))

names(pbs_mat)= c("chr", "mid", "AGA", "AVB", "CWB", "HAY", "TYD", "TYE", "WEA", "YHE")
#filter these for those that overlap genes

#bed file of trout genes from ensembl 
gene_ranges<-read_tsv("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/PicMin/gene_ranges.tsv")

pbs_mat$midPos<-as.numeric(pbs_mat$mid)

pbs_mat$start <- pbs_mat$midPos - 25000
pbs_mat$end <- pbs_mat$midPos + 25000

# transform into genomic ranges objects then overlap 
gr1<-makeGRangesFromDataFrame(gene_ranges,
                              keep.extra.columns = FALSE,
                              seqinfo = NULL,
                              seqnames.field=c("chr"),
                              start.field = "start", 
                              end.field = "end",
                              ignore.strand = TRUE,
                              na.rm=TRUE)

gr2<-makeGRangesFromDataFrame(pbs_mat,
                              keep.extra.columns = FALSE,
                              seqinfo = NULL,
                              seqnames.field=c("chr"),
                              start.field = "start", 
                              end.field = "end",
                              ignore.strand = TRUE,
                              na.rm=TRUE)

gr_sub<-subsetByOverlaps(gr2,gr1)
gr_dat_all_50k<-as.data.frame(gr_sub)
#save these windows
write_tsv(gr_dat_all_50k, "pbs_gene_windows_50k.tsv")

pbs_windows<-as.data.frame(read_tsv("pbs_matrix_genes_50k.tsv"))

pbs_windows$midPos <- (pbs_windows$start + 25000)

#create ID column then intersect pbs rows with these windows
pbs_mat$ID <- paste0(pbs_mat$chr, "_", pbs_mat$midPos)
pbs_windows$ID <- paste0(pbs_windows$seqnames, "_", pbs_windows$midPos)

pbs_res<-setDT(pbs_mat)[ID %chin%pbs_windows$ID]
pbs_genes<-as.data.frame(pbs_res)

#make the pbs values numeric - dunno why they weren't already 
pbs_genes$AGA = as.numeric(pbs_genes$AGA)
pbs_genes$AVB = as.numeric(pbs_genes$AVB)
pbs_genes$CWB = as.numeric(pbs_genes$CWB)
pbs_genes$HAY = as.numeric(pbs_genes$HAY)
pbs_genes$TYD = as.numeric(pbs_genes$TYD)
pbs_genes$TYE = as.numeric(pbs_genes$TYE)
pbs_genes$WEA = as.numeric(pbs_genes$WEA)
pbs_genes$YHE = as.numeric(pbs_genes$YHE)

#calc some bleddy empirical pvals (meaningless)
pbs_genes$AGA_pval<-qvalue::empPvals(stat = pbs_genes$AGA, stat0 = pbs_genes$AGA)
pbs_genes$AVB_pval<-qvalue::empPvals(stat = pbs_genes$AVB, stat0 = pbs_genes$AVB)
pbs_genes$CWB_pval<-qvalue::empPvals(stat = pbs_genes$CWB, stat0 = pbs_genes$CWB)
pbs_genes$HAY_pval<-qvalue::empPvals(stat = pbs_genes$HAY, stat0 = pbs_genes$HAY)
pbs_genes$TYD_pval<-qvalue::empPvals(stat = pbs_genes$TYD, stat0 = pbs_genes$TYD)
pbs_genes$TYE_pval<-qvalue::empPvals(stat = pbs_genes$TYE, stat0 = pbs_genes$TYE)
pbs_genes$WEA_pval<-qvalue::empPvals(stat = pbs_genes$WEA, stat0 = pbs_genes$WEA)
pbs_genes$YHE_pval<-qvalue::empPvals(stat = pbs_genes$YHE, stat0 = pbs_genes$YHE)

summary(pbs_genes$AGA_pval)

# okay lets extract the sites with longest branch stats
AGA_outliers<-pbs_genes[pbs_genes$AGA_pval<0.005,]
AGA_outliers <- subset(AGA_outliers, select = c("chr", "midPos", "start", "end", "AGA_pval"))

AVB_outliers<-pbs_genes[pbs_genes$AVB_pval<0.005,]
AVB_outliers <- subset(AVB_outliers, select = c("chr", "midPos", "start", "end", "AVB_pval"))

CWB_outliers<-pbs_genes[pbs_genes$CWB_pval<0.005,]
CWB_outliers <- subset(CWB_outliers, select = c("chr", "midPos", "start", "end", "CWB_pval"))

HAY_outliers<-pbs_genes[pbs_genes$HAY_pval<0.005,]
HAY_outliers <- subset(HAY_outliers, select = c("chr", "midPos", "start", "end", "HAY_pval"))

TYD_outliers<-pbs_genes[pbs_genes$TYD_pval<0.005,]
TYD_outliers <- subset(TYD_outliers, select = c("chr", "midPos", "start", "end", "TYD_pval"))

TYE_outliers<-pbs_genes[pbs_genes$TYE_pval<0.005,]
TYE_outliers <- subset(TYE_outliers, select = c("chr", "midPos", "start", "end", "TYE_pval"))

WEA_outliers<-pbs_genes[pbs_genes$WEA_pval<0.005,]
WEA_outliers <- subset(WEA_outliers, select = c("chr", "midPos", "start", "end", "WEA_pval"))

YHE_outliers<-pbs_genes[pbs_genes$YHE_pval<0.005,]
YHE_outliers <- subset(YHE_outliers, select = c("chr", "midPos", "start", "end", "YHE_pval"))

#ah bugger, this complicates the next step - let's have a rethink

outliers <- data.frame(subset(pbs_genes, AGA_pval < 0.005 | AVB_pval < 0.005 | CWB_pval < 0.005 | HAY_pval < 0.005 | TYD_pval < 0.005 | TYE_pval < 0.005 | WEA_pval < 0.005 | YHE_pval < 0.005))
outliers<- subset(outliers, select = c("chr", "midPos", "start", "end", "AGA_pval", "AVB_pval", "CWB_pval", "HAY_pval", "TYD_pval", "TYE_pval", "WEA_pval", "YHE_pval"))

# cool now let's see what genes sit at those outliers 

trout_genome = fread("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/PicMin/trout_genome_info.csv")
trout_genome = data.frame(subset(trout_genome, select = c("Name", "RefSeq")))
names(trout_genome) = c("name", "chr")
tophits = merge(outliers,trout_genome)

#window ids
tophits$window_id = paste0(tophits$name,":",as.integer(tophits$start),":",as.integer(tophits$end))

library("biomaRt")
# Now set up biomart
ensembl <- useMart("ensembl")
trout <- useDataset("strutta_gene_ensembl",mart=ensembl)

# Fetch genes
# got up to here but not grabbing genes - fix this # 
pbs_outlier_genes = getBM(attributes = c("ensembl_gene_id","chromosome_name","start_position","end_position","description","external_gene_name"),
                             filters = "chromosomal_region",
                             values = tophits$window_id,
                             mart=trout)
pbs_outlier_genes = data.table(pbs_outlier_genes)[order(chromosome_name,start_position),]
pbs_outlier_genes

gr1<- makeGRangesFromDataFrame(pbs_outlier_genes,
                               keep.extra.columns = TRUE,
                               seqinfo = NULL,
                               seqnames.field=c("chromosome_name"),
                               start.field = "start_position", 
                               end.field = "end_position",
                               ignore.strand = TRUE,
                               na.rm=TRUE)

gr2<- makeGRangesFromDataFrame(tophits,
                               keep.extra.columns = TRUE,
                               seqinfo = NULL,
                               seqnames.field=c("name"),
                               start.field = "start", 
                               end.field = "end",
                               ignore.strand = TRUE,
                               na.rm=TRUE)

library(data.table)
setkey(setDT(tophits), name, start, end)
setkey(setDT(pbs_outlier_genes), chromosome_name, start_position, end_position)
gene_annotated_pbs<-foverlaps(pbs_outlier_genes,tophits)
write_tsv(gene_annotated_pbs,"gene_annotated_pbs.tsv")
