library(ape)
library(MASS)
library(nlme)
library(lme4)
require(lmerTest)
library(optimx)
library(ggplot2)
library(MuMIn)
library(plotly)
library(plyr)
library(phytools)
library(scatterplot3d)
library(rgl)
library(gridExtra)

####read in pertinent data
#########
#df<-read.csv(file = "vertebrate_species_effects_no_GC_csv.txt", header = T)
#str(df)

df<-read.csv(file = "clean_ISD_effects.txt", header = T)
str(df)

#clean out unneccesary bits
df<-subset(df, select = -c(Error,df,P_value,t_value,df))
str(df)


dfGC<-read.csv(file = "CAAI_CAI_GC_1_29_2020.txt", header = T)
str(dfGC)

dfGC<-subset(dfGC, select = -c(Row,Speciesname))
str(dfGC)

#merge
df<-merge(df, dfGC, by.x="SpeciesName", by.y="NewickName")
str(df)

#read in order info
dfvert<-read.csv(file = "vertebrateage_order_sans_cow.txt", header = T)
str(dfvert)

#merge
df<-merge(df, dfvert, by.x="SpeciesName", by.y="Species_Name")
str(df)

#clean up
df<-subset(df, select = -c(ISD_correlation_coeff,CAI.y,CLUST_no_pfam_correlation_coeff,CLUST_pfam_correlation_coeff))
str(df)

########################set up subsets by order of life##########################

dfRod <-df[which(df$Order == "Rodentia"),names(df) %in% c("Species_Name","Estimate","fi_CAAI","Family","Order")]

dfCroc <-df[which(df$Order == "Crocodilia"),names(df) %in% c("Species_Name","Estimate","fi_CAAI","Family","Order")]

dfcow <-df[which(df$Order == "Artiodactyla"),names(df) %in% c("Species_Name","Estimate","fi_CAAI","Family","Order")]

dfPrim <-df[which(df$Order == "Primates"),names(df) %in% c("Species_Name","Estimate","fi_CAAI","Family","Order")]

dfCeta <-df[which(df$Order == "Cetacea"),names(df) %in% c("Species_Name","Estimate","fi_CAAI","Family","Order")]

dfCarn <-df[which(df$Order == "Carnivora"),names(df) %in% c("Species_Name","Estimate","fi_CAAI","Family","Order")]

dfchick <-df[which(df$Order == "Galliformes"),names(df) %in% c("Species_Name","Estimate","fi_CAAI","Family","Order")]

dfbat <-df[which(df$Order == "Chiroptera"),names(df) %in% c("Species_Name","Estimate","fi_CAAI","Family","Order")]

dfper <-df[which(df$Order == "Perciformes"),names(df) %in% c("Species_Name","Estimate","fi_CAAI","Family","Order")]

dfforplot<-rbind(dfRod, dfcow,dfPrim,dfbat)


##################################PLOT####################################

pvert <- ggplot(df, aes(fi_CAAI,Estimate)) + geom_point(color='black',size =2)+ labs( x = "CAIS", y = "Effect on ISD")  +geom_text(data = df, aes(x = 1.06, y = 0.018, label = "Spearman's R:0.66  \n Pearson's R:0.78\n p-value<2e-16"), colour = 'black', size = 7)
pvert <- pvert + stat_smooth(method = "lm", color="red")+theme(axis.title=element_text(size=28, face ="bold"),title=element_text(size=14, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))+ theme(axis.text = element_text(size = 24)) 

#pvert<-pvert + geom_point(data=dfCarn, aes(x=fi_CAAI,y=Estimate), color='green',size=3)
pvert<-pvert + geom_point(data=dfforplot, aes(x=fi_CAAI,y=Estimate, color= Order),size=4)+ theme(legend.position = c(0.8, 0.2))+theme(legend.text = element_text(colour="black", size=18, 
                                                                                                                                                                face="bold"))+ theme(legend.background = element_rect(fill="white",
                                                                                                                                                                                                                      size=0.5, linetype="solid", 
                                                                                                                                                                                                                      colour ="black"))


pvert

