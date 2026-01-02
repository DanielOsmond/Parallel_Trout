setwd("C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics/sfs/thetas")
library(readr)
library(ggplot2)#
library(dplyr)

TYD_thetas<-read.table("TYD.thetas.windows.gz.pestPG", header= F, sep = "\t")
head(TYD_thetas)
names(TYD_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                     "fuf", "fud", "fayh", "zeng", "nSites")

# ggplot(TYD_thetas_25) +
#   geom_jitter(aes(x=WinCenter, y=Tajima),
#             ) +
#   theme_bw()

TYE_thetas<-read.table("TYE.thetas.windows.gz.pestPG", header= F, sep = "\t")
names(TYE_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                     "fuf", "fud", "fayh", "zeng", "nSites")

TYA_thetas<-read.table("TYA.thetas.windows.gz.pestPG", header= F, sep = "\t")
names(TYA_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                     "fuf", "fud", "fayh", "zeng", "nSites")

HAY_thetas<-read.table("HAY1.thetas.windows.gz.pestPG", header= F, sep = "\t")
names(HAY_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                     "fuf", "fud", "fayh", "zeng", "nSites")
HAY_thetas$Pi <- HAY_thetas$tP / HAY_thetas$nSites

CAM_thetas<-read.table("CAM.thetas.windows.gz.pestPG", header= F, sep = "\t")
names(CAM_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                     "fuf", "fud", "fayh", "zeng", "nSites")

CWA_thetas<-read.table("CWA.thetas.windows.gz.pestPG", header= F, sep = "\t")
names(CWA_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                     "fuf", "fud", "fayh", "zeng", "nSites")
CWA_thetas$Pi <- CWA_thetas$tP / CWA_thetas$nSites

CWB_thetas<-read.table("CWB.thetas.windows.gz.pestPG", header= F, sep = "\t")
names(CWB_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                     "fuf", "fud", "fayh", "zeng", "nSites")
CWB_thetas$Pi <- CWB_thetas$tP / CWB_thetas$nSites

CWC_thetas<-read.table("CWC.thetas.windows.gz.pestPG", header= F, sep = "\t")
names(CWC_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                     "fuf", "fud", "fayh", "zeng", "nSites")
CWC_thetas$Pi <- CWC_thetas$tP / CWC_thetas$nSites


WEA_thetas<-read.table("WEA.thetas.windows.gz.pestPG", header= F, sep = "\t")
names(WEA_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                     "fuf", "fud", "fayh", "zeng", "nSites")
WEA_thetas$Pi <- WEA_thetas$tP / WEA_thetas$nSites

WEF_thetas<-read.table("WEF.thetas.windows.gz.pestPG", header= F, sep = "\t")
names(WEF_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                     "fuf", "fud", "fayh", "zeng", "nSites")
WEF_thetas$Pi <- WEF_thetas$tP / WEF_thetas$nSites

AGA_thetas<-read.table("AGA.thetas.windows.gz.pestPG", header= F, sep = "\t")
names(AGA_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                     "fuf", "fud", "fayh", "zeng", "nSites")
AGA_thetas$Pi <- AGA_thetas$tP / AGA_thetas$nSites

AGB_thetas<-read.table("AGB.thetas.windows.gz.pestPG", header= F, sep = "\t")
names(AGB_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                     "fuf", "fud", "fayh", "zeng", "nSites")
AGB_thetas$Pi <- AGB_thetas$tP / AGB_thetas$nSites

AVB_thetas<-read.table("AVB.thetas.windows.gz.pestPG", header= F, sep = "\t")
names(AVB_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                     "fuf", "fud", "fayh", "zeng", "nSites")
AVB_thetas$Pi <- AVB_thetas$tP / AVB_thetas$nSites

AVD_thetas<-read.table("AVD.thetas.windows.gz.pestPG", header= F, sep = "\t")
names(AVD_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                     "fuf", "fud", "fayh", "zeng", "nSites")
AVD_thetas$Pi <- AVD_thetas$tP / AVD_thetas$nSites

OWE_thetas<-read.table("OWE.thetas.windows.gz.pestPG", header= F, sep = "\t")
names(OWE_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                     "fuf", "fud", "fayh", "zeng", "nSites")
OWE_thetas$Pi <- OWE_thetas$tP / OWE_thetas$nSites

NEV2_thetas<-read.table("NEV2.thetas.windows.gz.pestPG", header= F, sep = "\t")
names(NEV2_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                     "fuf", "fud", "fayh", "zeng", "nSites")
NEV2_thetas$Pi <- NEV2_thetas$tP / NEV2_thetas$nSites

YHA_thetas<-read.table("YHA.thetas.windows.gz.pestPG", header= F, sep = "\t")
names(YHA_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                      "fuf", "fud", "fayh", "zeng", "nSites")
YHA_thetas$Pi <- YHA_thetas$tP / YHA_thetas$nSites

YHE_thetas<-read.table("YHE.thetas.windows.gz.pestPG", header= F, sep = "\t")
names(YHE_thetas)= c("metrics", "Chr", "WinCenter", "tW", "tP" ,"tF", "tH", "tL", "Tajima",
                      "fuf", "fud", "fayh", "zeng", "nSites")
YHE_thetas$Pi <- YHE_thetas$tP / YHE_thetas$nSites

CAM_thetas$Pi <- CAM_thetas$tP / CAM_thetas$nSites

TYD_thetas$Pi <- TYD_thetas$tP / TYD_thetas$nSites
TYE_thetas$Pi <- TYE_thetas$tP / TYE_thetas$nSites
TYA_thetas$Pi <- TYA_thetas$tP / TYA_thetas$nSites

Tyne_Tajima<-cbind.data.frame(TYD_thetas$Tajima, TYE_thetas$Tajima, TYA_thetas$Tajima)
Tyne_Pi<-cbind.data.frame(TYD_thetas$Pi, TYE_thetas$Pi, TYA_thetas$Pi)
names(Tyne_Tajima)=c("TYD", "TYE" ,"TYA")
names(Tyne_Pi)= c("TYD", "TYE" ,"TYA")
Tyne_Tajima_filtered <- subset(Tyne_Tajima, !TYD == 0 & !TYE == 0 & !TYA == 0)
Tyne_Pi_filtered <- subset(Tyne_Pi, !TYD == 0 & !TYE == 0 & !TYA == 0)

#### okay now find a way to plot all three densities on one graph
#

head(TYD_thetas)

TYA_thetas_25<- filter(TYA_thetas, Chr == "NC_042981.1") 
TYD_thetas_25<- filter(TYD_thetas, Chr == "NC_042981.1") 
TYE_thetas_25<- filter(TYE_thetas, Chr == "NC_042981.1") 
WEA_thetas_25<- filter(WEA_thetas, Chr == "NC_042981.1")
WEF_thetas_25<- filter(WEF_thetas, Chr == "NC_042981.1") 


TYE_taj_25<-ggplot(TYE_thetas_25)+
  geom_line(aes(x=WinCenter,y=Tajima))+
  scale_x_continuous(limits=c(20000000,30000000))+
  theme_bw()#

TYA_taj_25<-ggplot(TYA_thetas_25)+
  geom_line(aes(x=WinCenter,y=Tajima))+
  scale_x_continuous(limits=c(20000000,30000000))+
  theme_bw()

TYD_taj_25<-ggplot(TYD_thetas_25)+
  geom_line(aes(x=WinCenter,y=Tajima))+
  scale_x_continuous(limits=c(20000000,30000000))+
  theme_bw()

WEA_taj_25<-ggplot(WEA_thetas_25)+
  geom_line(aes(x=WinCenter,y=Tajima))+
  scale_x_continuous(limits=c(20000000,30000000))+
  theme_bw()

WEF_taj_25<-ggplot(WEF_thetas_25)+
  geom_line(aes(x=WinCenter,y=Tajima))+
  scale_x_continuous(limits=c(20000000,30000000))+
  theme_bw()

library(cowplot)
chr25_ne_taj<-plot_grid(TYA_taj_25, TYD_taj_25, TYE_taj_25, WEA_taj_25, WEF_taj_25,
          cols=1)
ggsave("chr25_ne_taj.png", width = 6, height =10)

TYE_Pi_25<-ggplot(TYE_thetas_25)+
  geom_line(aes(x=WinCenter,y=Pi))+
  scale_x_continuous(limits=c(20000000,30000000))+
  theme_bw()#

TYA_Pi_25<-ggplot(TYA_thetas_25)+
  geom_line(aes(x=WinCenter,y=Pi))+
  scale_x_continuous(limits=c(20000000,30000000))+
  theme_bw()

TYD_Pi_25<-ggplot(TYD_thetas_25)+
  geom_line(aes(x=WinCenter,y=Pi))+
  scale_x_continuous(limits=c(20000000,30000000))+
  theme_bw()

WEA_Pi_25<-ggplot(WEA_thetas_25)+
  geom_line(aes(x=WinCenter,y=Pi))+
  scale_x_continuous(limits=c(20000000,30000000))+
  theme_bw()

WEF_Pi_25<-ggplot(WEF_thetas_25)+
  geom_line(aes(x=WinCenter,y=Pi))+
  scale_x_continuous(limits=c(20000000,30000000))+
  theme_bw()

chr25_ne_Pi<-plot_grid(TYA_Pi_25, TYD_Pi_25, TYE_Pi_25, WEA_Pi_25, WEF_Pi_25,
                        cols=1)
ggsave("chr25_ne_Pi.png", width = 6, height =10)

ggplot(TYD_thetas)+
  geom_density(aes(x=Pi))+
  theme_bw()

ylim=c(0,0.7)
xlim=c(-2,3)
piy=c(0,5)
pix=c(0,0.5)

Tyne_Taj<-ggplot(data = Tyne_Tajima_filtered) +
  geom_density(aes(x = TYD, color = "TYD"), alpha = 0.5) +
  geom_density(aes(x = TYE, color = "TYE"), alpha = 0.5) +
  geom_density(aes(x = TYA, color = "TYA"), alpha = 0.5) +
  labs(x = "Tajima's D",
       y = "Density") +
  scale_color_manual(values = c("green", "orange", "deeppink")) +
  scale_y_continuous(limits = ylim)+
  scale_x_continuous(limits =xlim)+
  theme_bw() +
  theme(legend.position = "none")

#ggsave("Tyne_Tajima.png", width = 4, height =3)

Tyne_Pi<-ggplot(data = Tyne_Pi_filtered) +
  geom_density(aes(x = TYD, color = "TYD"), alpha = 0.5) +
  geom_density(aes(x = TYE, color = "TYE"), alpha = 0.5) +
  geom_density(aes(x = TYA, color = "TYA"), alpha = 0.5) +
  labs(x = "Pi",
       y = "Density") +
  scale_color_manual(values = c("green", "orange", "deeppink")) +
  scale_y_continuous(limits = piy)+
  scale_x_continuous(limits =pix)+
  theme_bw()
#ggsave("Tyne_Pi.png", width = 4, height =3)



Kernow_Tajima<-cbind.data.frame(HAY_thetas$Tajima, CAM_thetas$Tajima)
Kernow_Pi<-cbind.data.frame(HAY_thetas$Pi, CAM_thetas$Pi)
names(Kernow_Tajima)=c("HAY", "CAM")
names(Kernow_Pi)= c("HAY", "CAM")
Kernow_Tajima_filtered <- subset(Kernow_Tajima, !HAY == 0 & !CAM == 0)
Kernow_Pi_filtered <- subset(Kernow_Pi, !HAY == 0 & !CAM == 0)

Kernow_Taj<-ggplot(data = Kernow_Tajima_filtered) +
  geom_density(aes(x = HAY, color = "HAY"), alpha = 0.5) +
  geom_density(aes(x = CAM, color = "CAM"), alpha = 0.5) +
  labs(x = "Tajima's D",
       y = "Density") +
  scale_color_manual(values = c("green", "orange")) +
  scale_y_continuous(limits = ylim)+
  scale_x_continuous(limits =xlim)+
  theme_bw() +
  theme(legend.position = "none")

#ggsave("Kernow_Tajima.png", width = 4, height =3)

Kernow_Pi<-ggplot(data = Kernow_Pi_filtered) +
  geom_density(aes(x = HAY, color = "HAY"), alpha = 0.5) +
  geom_density(aes(x = CAM, color = "CAM"), alpha = 0.5) +
    labs(x = "Pi",
       y = "Density") +
  scale_color_manual(values = c("green", "orange")) +
  scale_y_continuous(limits = piy)+
  scale_x_continuous(limits =pix)+
  theme_bw()
#ggsave("Kernow_Pi.png", width = 4, height =3)



Conwy_Tajima<-cbind.data.frame(CWA_thetas$Tajima, CWB_thetas$Tajima, CWC_thetas$Tajima)
Conwy_Pi<-cbind.data.frame(CWA_thetas$Pi, CWB_thetas$Pi, CWC_thetas$Pi)
names(Conwy_Tajima)=c("CWA", "CWB", "CWC")
names(Conwy_Pi)= c("CWA", "CWB", "CWC")
Conwy_Tajima_filt<- subset(Conwy_Tajima, !CWA == 0 & !CWB == 0 & !CWC ==0)
Conwy_Pi_filtered <- subset(Conwy_Pi, !CWA == 0 & !CWB == 0 & !CWC == 0)

Conwy_Nev_Taj<-ggplot(data = Conwy_Tajima_filt) +
  geom_density(aes(x = CWA, color = "CWA"), alpha = 0.5) +
  geom_density(aes(x = CWB, color = "CWB"), alpha = 0.5) +
  geom_density(aes(x = CWC, color = "CWC"), alpha = 0.5) +
 # geom_density(aes(x = NEV, color = "NEV"), alpha = 0.5) +
  labs(x = "Tajima's D",
       y = "Density") +
 scale_color_manual(values = c("orange", "deeppink", "green")) +
  scale_y_continuous(limits = ylim)+
  scale_x_continuous(limits =xlim)+
  theme_bw() +
  theme(legend.position = "none")

#ggsave("Conwy_Nev_Tajima.png", width = 4, height =3)

Conwy_Nev_Pi<-ggplot(data = Conwy_Pi_filtered) +
  geom_density(aes(x = CWA, color = "CWA"), alpha = 0.5) +
  geom_density(aes(x = CWB, color = "CWB"), alpha = 0.5) +
  geom_density(aes(x = CWC, color = "CWC"), alpha = 0.5) +
  #geom_density(aes(x = NEV, color = "NEV"), alpha = 0.5) +
  labs(x = "Pi",
       y = "Density") +
  scale_color_manual(values = c("orange", "deeppink", "green")) +
  scale_y_continuous(limits = piy)+
  scale_x_continuous(limits =pix)+
  theme_bw()
#ggsave("Conwy_Nev_Pi.png", width = 4, height =3)


Ireland_Tajima<-cbind.data.frame(AVB_thetas$Tajima, AVD_thetas$Tajima, OWE_thetas$Tajima)
Ireland_Pi<-cbind.data.frame(AVB_thetas$Pi, AVD_thetas$Pi, OWE_thetas$Pi)
names(Ireland_Tajima)=c("AVB", "AVD", "OWE")
names(Ireland_Pi)= c("AVB", "AVD", "OWE")
Ireland_Tajima_filtered <- subset(Ireland_Tajima, !AVB == 0 & !AVD == 0 & !OWE == 0)
Ireland_Pi_filtered <- subset(Ireland_Pi, !AVB == 0 & !AVD == 0 & !OWE == 0)

Ireland_Taj<-ggplot(data = Ireland_Tajima_filtered) +
  geom_density(aes(x = AVB, color = "AVB"), alpha = 0.5) +
  geom_density(aes(x = AVD, color = "AVD"), alpha = 0.5) +
  geom_density(aes(x = OWE, color = "OWE"), alpha = 0.5) +
  labs(x = "Tajima's D",
       y = "Density") +
  scale_color_manual(values = c("orange", "deeppink", "green")) +
  scale_y_continuous(limits = ylim)+
  scale_x_continuous(limits =xlim)+
  theme_bw() +
  theme(legend.position = "none")
#ggsave("Ireland_Tajima.png", width = 4, height =3)

Ireland_Pi<-ggplot(data = Ireland_Pi_filtered) +
  geom_density(aes(x = AVB, color = "AVB"), alpha = 0.5) +
  geom_density(aes(x = AVD, color = "AVD"), alpha = 0.5) +
  geom_density(aes(x = OWE, color = "OWE"), alpha = 0.5) +
  labs(x = "Pi",
       y = "Density") +
  scale_color_manual(values = c("orange", "deeppink", "green")) +
  scale_y_continuous(limits = piy)+
  scale_x_continuous(limits =pix)+
  theme_bw()
#ggsave("Ireland_Pi.png", width = 4, height =3)

Anglesey_Tajima<-cbind.data.frame(AGA_thetas$Tajima, AGB_thetas$Tajima, NEV2_thetas$Tajima)
Anglesey_Pi<-cbind.data.frame(AGA_thetas$Pi, AGB_thetas$Pi, NEV2_thetas$Pi)
names(Anglesey_Tajima)=c("AGA", "AGB", "NEV")
names(Anglesey_Pi)= c("AGA", "AGB", "NEV")
Anglesey_Tajima_filtered <- subset(Anglesey_Tajima, !AGA == 0 & !AGB == 0 & !NEV == 0)
Anglesey_Pi_filtered <- subset(Anglesey_Pi, !AGA == 0 & !AGB == 0 & !NEV == 0)

Anglesey_Taj<-ggplot(data = Anglesey_Tajima_filtered) +
  geom_density(aes(x = AGA, color = "AGA"), alpha = 0.5) +
  geom_density(aes(x = AGB, color = "AGB"), alpha = 0.5) +
  geom_density(aes(x = NEV, color = "NEV"), alpha = 0.5) +
  labs(x = "Tajima's D",
       y = "Density") +
  scale_color_manual(values = c("orange", "deeppink", "green")) +
  scale_y_continuous(limits = ylim)+
  scale_x_continuous(limits =xlim)+
  theme_bw() +
  theme(legend.position = "none")
#ggsave("Anglesey_Tajima.png", width = 4, height =3)

Anglesey_Pi<-ggplot(data = Anglesey_Pi_filtered) +
  geom_density(aes(x = AGA, color = "AGA"), alpha = 0.5) +
  geom_density(aes(x = AGB, color = "AGB"), alpha = 0.5) +
  geom_density(aes(x = NEV, color = "NEV"), alpha = 0.5) +
  labs(x = "Pi",
       y = "Density") +
  scale_color_manual(values = c("orange", "deeppink", "green")) +
  scale_y_continuous(limits = piy)+
  scale_x_continuous(limits =pix)+
  theme_bw()
#ggsave("Anglesey_Pi.png", width = 4, height =3)

Wear_Tajima<-cbind.data.frame(WEA_thetas$Tajima, WEF_thetas$Tajima)
Wear_Pi<-cbind.data.frame(WEA_thetas$Pi, WEF_thetas$Pi)
names(Wear_Tajima)=c("WEA", "WEF")
names(Wear_Pi)= c("WEA", "WEF")
Wear_Tajima_filtered <- subset(Wear_Tajima, !WEA == 0 & !WEF == 0)
Wear_Pi_filtered <- subset(Wear_Pi, !WEA == 0 & !WEF == 0)

Wear_Taj<-ggplot(data = Wear_Tajima_filtered) +
  geom_density(aes(x = WEA, color = "WEA"), alpha = 0.5) +
  geom_density(aes(x = WEF, color = "WEF"), alpha = 0.5) +
  labs(x = "Tajima's D",
       y = "Density") +
  scale_color_manual(values = c("orange", "green")) +
  scale_y_continuous(limits = ylim)+
  scale_x_continuous(limits =xlim)+
  theme_bw() +
  theme(legend.position = "none")
Wear_Taj
#ggsave("Wear_Tajima.png", width = 4, height =3)

Wear_Pi<-ggplot(data = Wear_Pi_filtered) +
  geom_density(aes(x = WEA, color = "WEA"), alpha = 0.5) +
  geom_density(aes(x = WEF, color = "WEF"), alpha = 0.5) +
  labs(x = "Pi",
       y = "Density") +
  scale_color_manual(values = c("orange", "green")) +
  scale_y_continuous(limits = piy)+
  scale_x_continuous(limits =pix)+
  theme_bw()
Wear_Pi
#ggsave("Wear_Pi.png", width = 4, height =3)

Ystwyth_Tajima<-cbind.data.frame(YHA_thetas$Tajima, YHE_thetas$Tajima, NEV2_thetas$Tajima)
Ystwyth_Pi<-cbind.data.frame(YHA_thetas$Pi, YHE_thetas$Pi, NEV2_thetas$Pi)
names(Ystwyth_Tajima)=c("YHA", "YHE", "NEV")
names(Ystwyth_Pi)= c("YHA", "YHE", "NEV")
Ystwyth_Tajima_filtered <- subset(Ystwyth_Tajima, !YHA == 0 & !YHE == 0 & !NEV == 0)
Ystwyth_Pi_filtered <- subset(Ystwyth_Pi, !YHA == 0 & !YHE == 0 & !NEV == 0)

Ystwyth_Taj<-ggplot(data = Ystwyth_Tajima_filtered) +
  geom_density(aes(x = NEV, color = "NEV"), alpha = 0.5) +
  geom_density(aes(x = YHA, color = "YHA"), alpha = 0.5) +
  geom_density(aes(x = YHE, color = "YHE"), alpha = 0.5) +
  labs(x = "Tajima's D",
       y = "Density") +
  scale_color_manual(values = c("blue", "green", "orange")) +
  scale_y_continuous(limits = ylim)+
  scale_x_continuous(limits =xlim)+
  theme_bw() +
  theme(legend.position = "none") 
Ystwyth_Taj
#ggsave("Ystwyth_Tajima.png", width = 4, height =3)

Ystwyth_Pi<-ggplot(data = Ystwyth_Pi_filtered) +
  geom_density(aes(x = NEV, color = "NEV"), alpha = 0.5) +
  geom_density(aes(x = YHA, color = "YHA"), alpha = 0.5) +
  geom_density(aes(x = YHE, color = "YHE"), alpha = 0.5) +
  labs(x = "Pi",
       y = "Density") +
  scale_color_manual(values = c("blue", "green", "orange")) +
  scale_y_continuous(limits = piy)+
  scale_x_continuous(limits =pix)+
  theme_bw()
Ystwyth_Pi
#ggsave("Ystwyth_Pi.png", width = 4, height =3)


library(cowplot)
combined_nucs<- plot_grid(Anglesey_Taj, Anglesey_Pi, Ireland_Taj, Ireland_Pi, Kernow_Taj, Kernow_Pi,
          Tyne_Taj, Tyne_Pi, Wear_Taj, Wear_Pi, Conwy_Nev_Taj, Conwy_Nev_Pi, 
          Ystwyth_Taj, Ystwyth_Pi, ncol = 2, rel_widths = c(3,4))
ggsave("combined_pops_Taj_Pi.png", width = 8, height = 14)
