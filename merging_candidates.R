## An Rscript written by DROsmond 9.11.23 with a troubled mental state 
## Yes I know I should write functions - write these if I ever publicise this code 
## in the meantime, forgive a struggling end of PhD soul.

setwd("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics")

library(dplyr)
library(readr)
library(ggplot2)
library(tidyverse)
library(GenomicRanges)
library(biomaRt)

## read in the contigs to rename to chrs ##

contigs<-read.table("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/contigs_names.txt", header=F, sep="\t")
names(contigs) = c("name1", "name2")
name1<-as.vector(contigs$name1)
name2<-as.vector(contigs$name2)

## read in our files

Cu_eig1_99<- read.table("./mafs/Cu_eig1_99%.txt", header=T, sep = "\t")
Cu_eig1_99_fully<- read.table("./mafs/Cu_eig1_99%_fully_parallel.txt", header=T, sep = "\t")

cymru_eig1_99<-read.table("./mafs/cymru_eig1_99%.txt", header = T, sep = "\t")
cymru_eig1_99_fully<-read.table("./mafs/cymru_eig1_99%_fully_parallel.txt", header = T, sep = "\t")

NE_eig1_99<- read.table("./mafs/NE_eig1_99%.txt", header = T, sep = "\t")
NE_eig1_99_fully<-read.table("./mafs/NE_eig1_99%_fully_parallel.txt", header = T, sep = "\t")

Pb_eig1_99 <- read.table("./mafs/Pb_eig1_99%.txt", header = T, sep = "\t")
Pb_eig1_99_fully <- read.table("./mafs/Pb_eig1_99%_fully_parallel.txt", header = T, sep = "\t")

UK_eig1_99 <- read.table("./mafs/UK_eig1_99%.txt", header = T, sep = "\t")
UK_eig1_99_fully <- read.table("./mafs/UK_eig1_99%_fully_parallel.txt", header = T, sep = "\t")

### now fix the dataframes for overlap

Cu_eig1_99 <- separate(Cu_eig1_99, window_id, into = c("chr", "pos"), sep = ":", remove = T)
Cu_eig1_99 <- separate(Cu_eig1_99, pos, into = c("start", "end"), sep = "-", remove = T)
Cu_eig1_99 <- merge(Cu_eig1_99, contigs, by.x = "chr", by.y = "name1", all.x = TRUE)
Cu_eig1_99$chr <- Cu_eig1_99$name2
Cu_eig1_99$chr<-as.numeric(Cu_eig1_99$chr)

Pb_eig1_99 <- separate(Pb_eig1_99, window_id, into = c("chr", "pos"), sep = ":", remove = T)
Pb_eig1_99 <- separate(Pb_eig1_99, pos, into = c("start", "end"), sep = "-", remove = T)
Pb_eig1_99 <- merge(Pb_eig1_99, contigs, by.x = "chr", by.y = "name1", all.x = TRUE)
Pb_eig1_99$chr <- Pb_eig1_99$name2
Pb_eig1_99$chr<-as.numeric(Pb_eig1_99$chr)

cymru_eig1_99 <- separate(cymru_eig1_99, window_id, into = c("chr", "pos"), sep = ":", remove = T)
cymru_eig1_99 <- separate(cymru_eig1_99, pos, into = c("start", "end"), sep = "-", remove = T)
cymru_eig1_99$start <- as.numeric(cymru_eig1_99$start)
cymru_eig1_99$end <- as.numeric(cymru_eig1_99$end)
cymru_eig1_99 <- merge(cymru_eig1_99, contigs, by.x = "chr", by.y = "name1", all.x = TRUE)
cymru_eig1_99$chr <- cymru_eig1_99$name2
cymru_eig1_99$chr<-as.numeric(cymru_eig1_99$chr)

NE_eig1_99 <- separate(NE_eig1_99, window_id, into = c("chr", "pos"), sep = ":", remove = T)
NE_eig1_99 <- separate(NE_eig1_99, pos, into = c("start", "end"), sep = "-", remove = T)
NE_eig1_99 <- merge(NE_eig1_99, contigs, by.x = "chr", by.y = "name1", all.x = TRUE)
NE_eig1_99$chr <- NE_eig1_99$name2
NE_eig1_99$chr<-as.numeric(NE_eig1_99$chr)

UK_eig1_99 <- separate(UK_eig1_99, window_id, into = c("chr", "pos"), sep = ":", remove = T)
UK_eig1_99 <- separate(UK_eig1_99, pos, into = c("start", "end"), sep = "-", remove = T)
UK_eig1_99 <- merge(UK_eig1_99, contigs, by.x = "chr", by.y = "name1", all.x = TRUE)
UK_eig1_99$chr <- UK_eig1_99$name2
UK_eig1_99$chr<-as.numeric(UK_eig1_99$chr)

Cu_eig1_99_fully <- separate(Cu_eig1_99_fully, window_id, into = c("chr", "pos"), sep = ":", remove = T)
Cu_eig1_99_fully <- separate(Cu_eig1_99_fully, pos, into = c("start", "end"), sep = "-", remove = T)
Cu_eig1_99_fully <- merge(Cu_eig1_99_fully, contigs, by.x = "chr", by.y = "name1", all.x = TRUE)
Cu_eig1_99_fully$chr <- Cu_eig1_99_fully$name2
Cu_eig1_99_fully$chr<-as.numeric(Cu_eig1_99_fully$chr)

Pb_eig1_99_fully <- separate(Pb_eig1_99_fully, window_id, into = c("chr", "pos"), sep = ":", remove = T)
Pb_eig1_99_fully <- separate(Pb_eig1_99_fully, pos, into = c("start", "end"), sep = "-", remove = T)
Pb_eig1_99_fully <- merge(Pb_eig1_99_fully, contigs, by.x = "chr", by.y = "name1", all.x = TRUE)
Pb_eig1_99_fully$chr <- Pb_eig1_99_fully$name2
Pb_eig1_99_fully$chr<-as.numeric(Pb_eig1_99_fully$chr)

cymru_eig1_99_fully <- separate(cymru_eig1_99_fully, window_id, into = c("chr", "pos"), sep = ":", remove = T)
cymru_eig1_99_fully <- separate(cymru_eig1_99_fully, pos, into = c("start", "end"), sep = "-", remove = T)
cymru_eig1_99_fully <- merge(cymru_eig1_99_fully, contigs, by.x = "chr", by.y = "name1", all.x = TRUE)
cymru_eig1_99_fully$chr <- cymru_eig1_99_fully$name2
cymru_eig1_99_fully$chr<-as.numeric(cymru_eig1_99_fully$chr)

NE_eig1_99_fully <- separate(NE_eig1_99_fully, window_id, into = c("chr", "pos"), sep = ":", remove = T)
NE_eig1_99_fully <- separate(NE_eig1_99_fully, pos, into = c("start", "end"), sep = "-", remove = T)
NE_eig1_99_fully <- merge(NE_eig1_99_fully, contigs, by.x = "chr", by.y = "name1", all.x = TRUE)
NE_eig1_99_fully$chr <- NE_eig1_99_fully$name2
NE_eig1_99_fully$chr<-as.numeric(NE_eig1_99_fully$chr)

UK_eig1_99_fully <- separate(UK_eig1_99_fully, window_id, into = c("chr", "pos"), sep = ":", remove = T)
UK_eig1_99_fully <- separate(UK_eig1_99_fully, pos, into = c("start", "end"), sep = "-", remove = T)
UK_eig1_99_fully <- merge(UK_eig1_99_fully, contigs, by.x = "chr", by.y = "name1", all.x = TRUE)
UK_eig1_99_fully$chr <- UK_eig1_99_fully$name2
UK_eig1_99_fully$chr<-as.numeric(UK_eig1_99_fully$chr)

### okay now read in the picmin table ### 

picmin05<-read.table("./PicMin/picmin_fdr05_contribs.tsv", header = T, sep = "\t")
picmin05$start<-picmin05$mid -25000
picmin05$end<-picmin05$mid+25000

picmin05 <- merge(picmin05, contigs, by.x = "chr", by.y = "name1", all.x = TRUE)
picmin05$chr <- picmin05$name2rele
picmin05$chr<-as.numeric(picmin05$chr)

## make some overlaps by genomic ranges ## 
## UKwide AF_vapeR set

picmin05_gr<-makeGRangesFromDataFrame(picmin05,
                                 keep.extra.columns = TRUE,rm 
                                 seqinfo = NULL,
                                 seqnames.field=c("chr"),
                                 start.field = "start", 
                                 end.field = "end",
                                 ignore.strand = TRUE,
                                 na.rm=TRUE)

UK_eig1_gr<-makeGRangesFromDataFrame(UK_eig1_99,
                                 keep.extra.columns = TRUE,
                                 seqinfo = NULL,
                                 seqnames.field=c("chr"),
                                 start.field = "start", 
                                 end.field = "end",
                                 ignore.strand = TRUE,
                                 na.rm=TRUE)

UK_eig1_picmin05<-subsetByOverlaps(picmin05_gr, UK_eig1_gr)
UK.AF_picmin05_overlap <-as.data.frame(UK_eig1_picmin05)
write_tsv(UK.AF_picmin05_overlap, "UK.AF_picmin05_overlap.tsv")

UK_eig1_gr<-makeGRangesFromDataFrame(UK_eig1_99_fully,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

UK_eig1_picmin05<-subsetByOverlaps(picmin05_gr, UK_eig1_gr)
UK.AF.fully_picmin05_overlap <-as.data.frame(UK_eig1_picmin05)
write_tsv(UK.AF.fully_picmin05_overlap, "UK.AF.fully_picmin05_overlap.tsv")

### NE AF-vapeR ### 

NE_eig1_gr<-makeGRangesFromDataFrame(NE_eig1_99,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

NE_eig1_picmin05<-subsetByOverlaps(picmin05_gr, NE_eig1_gr)
NE.AF_picmin05_overlap <-as.data.frame(NE_eig1_picmin05)
NE.AF_picmin05_overlap <- subset(NE.AF_picmin05_overlap, Tyne1_p <=0.05 & Wear_p <=0.05)
write_tsv(NE.AF_picmin05_overlap, "NE.AF_picmin05_overlap.tsv")

NE_eig1_gr<-makeGRangesFromDataFrame(NE_eig1_99_fully,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

NE_eig1_picmin05<-subsetByOverlaps(picmin05_gr, NE_eig1_gr)
NE.AF.fully_picmin05_overlap <-as.data.frame(NE_eig1_picmin05)
NE.AF.fully_picmin05_overlap <- subset(NE.AF.fully_picmin05_overlap, Tyne1_p <=0.05 & Wear_p <=0.05)
write_tsv(NE.AF.fully_picmin05_overlap, "NE.AF.fully_picmin05_overlap.tsv")

### Cu AF-vapeR ### 

Cu_eig1_gr<-makeGRangesFromDataFrame(Cu_eig1_99,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

Cu_eig1_picmin05<-subsetByOverlaps(picmin05_gr, Cu_eig1_gr)
Cu.AF_picmin05_overlap <-as.data.frame(Cu_eig1_picmin05)
Cu.AF_picmin05_overlap <- subset(Cu.AF_picmin05_overlap, Kernow_p <=0.05 & Anglesey_p <=0.05 & Ireland2_p <=0.05)
write_tsv(Cu.AF_picmin05_overlap, "Cu.AF_picmin05_overlap.tsv")

Cu_eig1_gr<-makeGRangesFromDataFrame(Cu_eig1_99_fully,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

Cu_eig1_picmin05<-subsetByOverlaps(picmin05_gr, Cu_eig1_gr)
Cu.AF.fully_picmin05_overlap <-as.data.frame(Cu_eig1_picmin05)
Cu.AF.fully_picmin05_overlap <- subset(Cu.AF.fully_picmin05_overlap, Kernow_p <=0.05 & Anglesey_p <=0.05 & Ireland2_p <=0.05)
write_tsv(Cu.AF.fully_picmin05_overlap, "Cu.AF.fully_picmin05_overlap.tsv")

## ### Pb AF-vapeR ### 

Pb_eig1_gr<-makeGRangesFromDataFrame(Pb_eig1_99,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

Pb_eig1_picmin05<-subsetByOverlaps(picmin05_gr, Pb_eig1_gr)
Pb.AF_picmin05_overlap <-as.data.frame(Pb_eig1_picmin05)
Pb.AF_picmin05_overlap <- subset(Pb.AF_picmin05_overlap, Nwales1_p <=0.05 & Tyne1_p <=0.05 & Wear_p <=0.05 & MidWales_p<=0.05 )
write_tsv(Pb.AF_picmin05_overlap, "Pb.AF_picmin05_overlap.tsv")

Pb_eig1_gr<-makeGRangesFromDataFrame(Pb_eig1_99_fully,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

Pb_eig1_picmin05<-subsetByOverlaps(picmin05_gr, Pb_eig1_gr)
Pb.AF.fully_picmin05_overlap <-as.data.frame(Pb_eig1_picmin05)
Pb.AF.fully_picmin05_overlap <- subset(Pb.AF.fully_picmin05_overlap, Nwales1_p <=0.05 & Tyne1_p <=0.05 & Wear_p <=0.05 & MidWales_p<=0.05 )
write_tsv(Pb.AF.fully_picmin05_overlap, "Pb.AF.fully_picmin05_overlap.tsv")

## ### Cymru AF-vapeR ### 

cymru_eig1_gr<-makeGRangesFromDataFrame(cymru_eig1_99,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

cymru_eig1_picmin05<-subsetByOverlaps(picmin05_gr, cymru_eig1_gr)
cymru.AF_picmin05_overlap <-as.data.frame(cymru_eig1_picmin05)
cymru.AF_picmin05_overlap <- subset(cymru.AF_picmin05_overlap, Nwales1_p <=0.05 & Anglesey_p <=0.05 & MidWales_p <=0.05)
write_tsv(cymru.AF_picmin05_overlap, "cymru.AF_picmin05_overlap.tsv")

cymru_eig1_gr<-makeGRangesFromDataFrame(cymru_eig1_99_fully,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

cymru_eig1_picmin05<-subsetByOverlaps(picmin05_gr, cymru_eig1_gr)
cymru.AF.fully_picmin05_overlap <-as.data.frame(cymru_eig1_picmin05)
cymru.AF.fully_picmin05_overlap <- subset(cymru.AF.fully_picmin05_overlap, Nwales1_p <=0.05 & Anglesey_p <=0.05 & MidWales_p <=0.05)
write_tsv(cymru.AF.fully_picmin05_overlap, "cymru.AF.fully_picmin05_overlap.tsv")

### pbs overlap ### 
pbs<-read_tsv("./sfs/pbs/gene_annotated_pbs.txt")

pbs_gr<-makeGRangesFromDataFrame(pbs,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chromosome_name"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

UK_eig1_gr<-makeGRangesFromDataFrame(UK_eig1_99,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

### overlap with UK AF_vapeR

UK_eig1_pbs<-subsetByOverlaps(pbs_gr, UK_eig1_gr)
UK.AF_pbs_overlap <-as.data.frame(UK_eig1_pbs)

## okay now just fully parallel

UK_eig1_gr<-makeGRangesFromDataFrame(UK_eig1_99_fully,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

UK_eig1_fully_pbs<-subsetByOverlaps(pbs_gr, UK_eig1_gr)
UK.AF.fully_pbs_overlap <-as.data.frame(UK_eig1_fully_pbs)
write_tsv(UK.AF.fully_pbs_overlap, "UK.AF.fully_pbs_overlap.tsv")

## copper ## 

Cu_eig1_gr<-makeGRangesFromDataFrame(Cu_eig1_99,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

Cu_eig1_pbs<-subsetByOverlaps(pbs_gr, Cu_eig1_gr)
Cu.AF_pbs_overlap <-as.data.frame(Cu_eig1_pbs)
Cu.AF_pbs_overlap <- subset(Cu.AF_pbs_overlap, AVB_pval <=0.05 | HAY_pval <=0.05 | AGA_pval <=0.05)
write_tsv(Cu.AF_pbs_overlap, "Cu.AF_pbs_overlap.tsv")

Cu_eig1_gr<-makeGRangesFromDataFrame(Cu_eig1_99_fully,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

Cu_eig1_pbs<-subsetByOverlaps(pbs_gr, Cu_eig1_gr)
Cu.AF_fully_pbs_overlap <-as.data.frame(Cu_eig1_pbs)
Cu.AF_pbs_overlap <- subset(Cu.AF_fully_pbs_overlap, AVB_pval <=0.05 | HAY_pval <=0.05 | AGA_pval <=0.05)
write_tsv(Cu.AF_pbs_overlap, "Cu.AF_pbs_overlap.tsv")

### lead ### 

Pb_eig1_gr<-makeGRangesFromDataFrame(Pb_eig1_99,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

Pb_eig1_pbs<-subsetByOverlaps(pbs_gr, Pb_eig1_gr)
Pb.AF_pbs_overlap <-as.data.frame(Pb_eig1_pbs)
Pb.AF_pbs_overlap <- subset(Pb.AF_pbs_overlap, CWB_pval <=0.05 | YHE_pval <=0.05 | TYD_pval <=0.05 | WEA_pval <= 0.05)
write_tsv(Pb.AF_pbs_overlap, "Pb.AF_pbs_overlap.tsv")

Pb_eig1_gr<-makeGRangesFromDataFrame(Pb_eig1_99_fully,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

Pb_eig1_pbs<-subsetByOverlaps(pbs_gr, Pb_eig1_gr)
Pb.AF_fully_pbs_overlap <-as.data.frame(Pb_eig1_pbs)
Pb.AF_pbs_overlap <- subset(Pb.AF_fully_pbs_overlap, CWB_pval <=0.05 | YHE_pval <=0.05 | TYD_pval <=0.05 | WEA_pval <= 0.05)
write_tsv(Pb.AF_pbs_overlap, "Pb.AF_fully_pbs_overlap.tsv")


### NE ### 

NE_eig1_gr<-makeGRangesFromDataFrame(NE_eig1_99,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

NE_eig1_pbs<-subsetByOverlaps(pbs_gr, NE_eig1_gr)
NE.AF_pbs_overlap <-as.data.frame(NE_eig1_pbs)
NE.AF_pbs_overlap <- subset(NE.AF_pbs_overlap, TYD_pval <=0.05 | TYE_pval <=0.05 | WEA_pval <=0.05)
write_tsv(NE.AF_pbs_overlap, "NE.AF_pbs_overlap.tsv")

NE_eig1_gr<-makeGRangesFromDataFrame(NE_eig1_99_fully,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

NE_eig1_pbs<-subsetByOverlaps(pbs_gr, NE_eig1_gr)
NE.AF_fully_pbs_overlap <-as.data.frame(NE_eig1_pbs)
NE.AF_pbs_overlap <- subset(NE.AF_fully_pbs_overlap, TYD_pval <=0.05 | TYE_pval <=0.05 | WEA_pval <=0.05)
write_tsv(NE.AF_pbs_overlap, "NE.AF_fully_pbs_overlap.tsv")

### NE ### 

cymru_eig1_gr<-makeGRangesFromDataFrame(cymru_eig1_99,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

cymru_eig1_pbs<-subsetByOverlaps(pbs_gr, cymru_eig1_gr)
cymru.AF_pbs_overlap <-as.data.frame(cymru_eig1_pbs)
cymru.AF_pbs_overlap <- subset(cymru.AF_pbs_overlap, YHE_pval <=0.05 | CWB_pval <=0.05 | AGA_pval <=0.05)
write_tsv(cymru.AF_pbs_overlap, "cymru.AF_pbs_overlap.tsv")

cymru_eig1_gr<-makeGRangesFromDataFrame(cymru_eig1_99_fully,
                                     keep.extra.columns = TRUE,
                                     seqinfo = NULL,
                                     seqnames.field=c("chr"),
                                     start.field = "start", 
                                     end.field = "end",
                                     ignore.strand = TRUE,
                                     na.rm=TRUE)

cymru_eig1_pbs<-subsetByOverlaps(pbs_gr, cymru_eig1_gr)
cymru.AF_fully_pbs_overlap <-as.data.frame(cymru_eig1_pbs)
cymru.AF_pbs_overlap <- subset(cymru.AF_fully_pbs_overlap, YHE_pval <=0.05 | CWB_pval <=0.05 | AGA_pval <=0.05)
write_tsv(cymru.AF_pbs_overlap, "cymru.AF_fully_pbs_overlap.tsv")
