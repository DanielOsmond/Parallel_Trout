#install.packages("UpSetR")
#install.packages("readr")
library("UpSetR")
library("readr")
library("magrittr")
library("dplyr")
library("tidyverse")

setwd("D:/phd/upset")

AF<-read_tsv("AF_vector.tsv")

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

AF$site <- paste0(AF$chr,"_",AF$pos)
Anglesey_poly <- subset(AF, AGA>0.025 & AGA<0.975 | AGB>0.025 & AGB<0.975| AGA<0.025 & AGB>0.975 | AGB<0.025 & AGA>0.975) ## retain polymophic sites, use 0.025 and 0.975 due to GL. 
Irish_1_poly  <- subset(AF, AVB>0.025 & AVB<0.975 | OWE>0.025 & OWE<0.975| AVB<0.025 & OWE>0.975 | OWE<0.025 & AVB>0.975) ## retain polymophic sites, use 0.025 and 0.975 due to GL.
Irish_2_poly  <- subset(AF, AVB>0.025 & AVB<0.975 | AVD>0.025 & AVD<0.975| AVD<0.025 & AVB>0.975 | AVB<0.025 & AVD>0.975)
Hayle_poly    <- subset(AF, HAY1>0.025 & HAY1<0.975 | CAM>0.025 & CAM<0.975| HAY1<0.025 & CAM>0.975 | CAM<0.025 & HAY1>0.975)
Mid_Wales_poly<- subset(AF, YHE>0.025 & YHE<0.975 | NEV2>0.025 & NEV2<0.975| YHE<0.025 & NEV2>0.975 | NEV2<0.025 & YHE>0.975)
N_Wales_poly  <- subset(AF, CWA>0.025 & CWA<0.975 | CWC>0.025 & CWC<0.975| CWA<0.025 & CWC>0.975 | CWC<0.025 & CWA>0.975)
N_Wales_2_poly<- subset(AF, CWB>0.025 & CWB<0.975 | CWC>0.025 & CWC<0.975| CWB<0.025 & CWC>0.975 | CWC<0.025 & CWB>0.975)
Tyne_1_poly   <- subset(AF, TYA>0.025 & TYA<0.975 | TYD>0.025 & TYD<0.975| TYA<0.025 & TYD>0.975 | TYD<0.025 & TYA>0.975)  
Tyne_2_poly   <- subset(AF, TYA>0.025 & TYA<0.975 | TYE>0.025 & TYE<0.975| TYA<0.025 & TYE>0.975 | TYE<0.025 & TYA>0.975)
Wear_poly     <- subset(AF, WEA>0.025 & WEA<0.975 | WEF>0.025 & WEF<0.975| WEA<0.025 & WEF>0.975 | WEF<0.025 & WEA>0.975)



summary(Anglesey_poly$AGA[Anglesey_poly$AGB==0])
summary(Anglesey_poly$AGA)

list_input <- list(Anglesey_poly$site, Irish_1_poly$site, Irish_2_poly$site, Hayle_poly$site,
                   Mid_Wales_poly$site, N_Wales_poly$site, N_Wales_2_poly$site,
                   Tyne_1_poly$site, Tyne_2_poly$site, Wear_poly$site) # just three pops to test function

UK_list <- list(Anglesey_poly$site, Irish_1_poly$site, Hayle_poly$site,
                Mid_Wales_poly$site, N_Wales_poly$site,
                Tyne_1_poly$site, Wear_poly$site)
names(list_input) <- c("Anglesey", "Irish1", "Irish2", "Cornwall", "Mid-Wales", "North-Wales1",
                       "North-Wales2", "Tyne1", "Tyne2", "Wear")
names(UK_list) <- c("Anglesey", "Irish",  "Cornwall", "Mid-Wales", "North-Wales",
                    "Tyne1", "Wear")
upsetplot<-upset(fromList(UK_list), order.by="freq", nsets=7, number.angles = 90,
                 mainbar.y.label = "Shared loci", sets.x.label = "polymorphic loci",
                 text.scale = c(1.3, 1.3, 1, 1, 2, 0.75))

pdf(file="UK_upset_freq.pdf")
upsetplot
dev.off()
AGAAGA

AF$site <- paste0(AF$chr,"_",AF$pos)
AGA_poly<-subset(AF, AGA>0.025 & AGA <0.975)
AGB_poly<-subset(AF, AGB>0.025 & AGB <0.975)
AVB_poly<-subset(AF, AVB>0.025 & AVB <0.975)
AVD_poly<-subset(AF, AVD>0.025 & AVD <0.975)
OWE_poly<-subset(AF, OWE>0.025 & OWE <0.975)
HAY1_poly<-subset(AF, HAY1>0.025 & HAY1 <0.975)
CAM_poly<-subset(AF, CAM>0.025 & CAM <0.975)
YHE_poly<-subset(AF, YHE>0.025 & YHE <0.975)
NEV2_poly<-subset(AF, NEV2>0.025 & NEV2 <0.975)
CWA_poly<-subset(AF, CWA>0.025 & CWA <0.975)
CWB_poly<-subset(AF, CWB>0.025 & CWB <0.975)
CWC_poly<-subset(AF, CWC>0.025 & CWC <0.975)
TYA_poly<-subset(AF, TYA>0.025 & TYA <0.975)
TYD_poly<-subset(AF, TYD>0.025 & TYD <0.975)
TYE_poly<-subset(AF, TYE>0.025 & TYE <0.975)
WEA_poly<-subset(AF, WEA>0.025 & WEA <0.975)
WEF_poly<-subset(AF, WEF>0.025 & WEF <0.975)

list_input <- list(AGA_poly$site, AGB_poly$site, AVB_poly$site, AVD_poly$site, OWE_poly$site,
                   HAY1_poly$site, CAM_poly$site,YHE_poly$site, NEV2_poly$site, CWA_poly$site,
                   CWB_poly$site, CWC_poly$site, TYA_poly$site, TYD_poly$site, TYE_poly$site,
                   WEA_poly$site, WEF_poly$site)
names(list_input) <- c("AGA", "AGB", "AVB", "AVD", "OWE", "HAY1", "CAM", "YHE", "NEV2", "CWA",
                       "CWB", "CWC", "TYA", "TYD", "TYE", "WEA", "WEF")

upsetplot<-upset(fromList(list_input), order.by="freq", 
                 sets=c("AGA", "AGB", "AVB", "AVD", "OWE", "HAY1", "CAM", "YHE", "NEV2", "CWA",
                        "CWB", "CWC", "TYA", "TYD", "TYE", "WEA", "WEF"), nsets=1, number.angles = 90,
                 mainbar.y.label = "Shared loci", sets.x.label = "polymorphic loci",
                 text.scale = c(1.3, 1.3, 1, 1, 2, 0.75))

pdf(file="UK_per_pop_upset_freq.pdf")
upsetplot
dev.off()

?upsetplot
