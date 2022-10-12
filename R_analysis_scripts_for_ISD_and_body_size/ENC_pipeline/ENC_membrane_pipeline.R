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

#read in newick tree
tree<-read.tree(file = "PhylogeneticTree_AllSpecies.nwk.txt")


####read in pertinent data
##########
dfvert<-read.csv(file = "clean_ISD_effects.txt", header = T)
str(dfvert)

dfverttrans<-read.csv(file = "clean_transmembrane_ISD_effects.txt", header = T)
str(dfverttrans)

dfvertnontrans<-read.csv(file = "clean_nontransmembrane_ISD_effects.txt", header = T)
str(dfvertnontrans)

#clean out unneccesary bits
dfvert<-subset(dfvert, select = -c(df,t_value))
str(dfvert)

dfGC<-read.csv(file = "CAAI_CAI_GC_1_29_2020.txt", header = T)
str(dfGC)

dfGC<-subset(dfGC, select = -c(Row,Speciesname))
str(dfGC)

dfENC<-read.csv(file = "ENC_6_8_2020.txt",header = T)
str(dfENC)

#merge
dfGC<-merge(dfGC, dfENC, by.x="SpeciesUID", by.y="SpeciesUID")
str(dfGC)

df<-merge(dfvert, dfGC, by.x="SpeciesName", by.y="NewickName")
str(df)

df<-merge(dfverttrans, df, by.x="SpeciesName", by.y="SpeciesName")
str(df)

df<-merge(dfvertnontrans, df, by.x="SpeciesName", by.y="SpeciesName")
str(df)

#########Comaparing dataset effects##########

lmtest1<-lm(df$Estimate~df$nontrans_Estimate)
summary(lmtest1)
cor.test(df$Estimate,df$nontrans_Estimate, method = 'pearson')

ISD1 <- plot(df$Estimate,df$nontrans_Estimate, xlab= "total ISD species effect", ylab="nontransmembrane ISD species effect") 

lmtest2<-lm(df$Estimate~df$trans_Estimate)
summary(lmtest2)
cor.test(df$Estimate,df$trans_Estimate, method = 'pearson')

ISD2 <- plot(df$Estimate,df$trans_Estimate, xlab= "total ISD species effect", ylab="transmembrane ISD species effect") 

lmtest3<-lm(df$trans_Estimate~df$nontrans_Estimate)
summary(lmtest3)
cor.test(df$trans_Estimate,df$nontrans_Estimate, method = 'spearman')

ISD3 <- ggplot(data = df, aes(trans_Estimate,nontrans_Estimate)) + labs(x= "Effect on ISD \n Transmembrane Pfams", y="Effect on ISD \n Nontransmembrane Pfams") + 
  geom_point(size=2)+ stat_smooth(method = "lm", color="red")+theme(axis.title=element_text(size=24, face ="bold"),title=element_text(size=20, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))+
  geom_text(data = df, aes(x = -0.005, y =0.018, label = " Spearman's R:0.59 \n Pearson's R:0.39 \n p-value=2.47e-12"),size=6, colour = 'black')+ theme(axis.text = element_text(size = 20)) +
  geom_errorbar(aes(ymin=nontrans_Estimate-nontrans_Std._Error, ymax=nontrans_Estimate+nontrans_Std._Error), colour="black")+coord_fixed(ratio=1)+
  geom_errorbarh(aes(xmin=trans_Estimate-trans_Std._Error, xmax=trans_Estimate+trans_Std._Error), colour="black") +
  annotate(geom = "segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
  annotate(geom = "segment", x = -Inf, xend = -Inf, y = -Inf, yend = -0.032,
           linetype = "dashed", color = "white")+
  annotate(geom = "segment", x = -Inf, xend = -Inf, y =  -0.03, yend = -0.01,
           linetype = "dashed", color = "white")+scale_y_continuous(limits=c(-0.01,0.02)) 
ISD3

membrane_error_plot<-ggplot()+geom_histogram(data=df, aes(trans_Std._Error, fill="Transmembrane"),binwidth = 0.00002,color = "black")+
  geom_histogram(data=df, aes(nontrans_Std._Error,fill="Nontransmembrane"),binwidth = 0.00002,color = "black")+
  labs(x= "Standard Error of Species Effect on ISD", y= "Count")+xlim(0,0.0018)+
  theme(legend.justification=c(1,0), legend.position=c(1,0.78))+theme(legend.title = element_text(colour="black", size=18, face="bold"),legend.text = element_text(colour="black", size=16, face="bold"))+
  labs(fill='Subset')+
  theme(axis.title=element_text(size=24, face ="bold"),axis.text=element_text(size=20, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
membrane_error_plot

mean(df$nontrans_Std._Error)
######################RAW#############

lmish<-lm(df$trans_Estimate~-df$Nc_nov, weights=(1/(df$trans_Std._Error))^2)
summary(lmish)
lmish<-lm(df$trans_Estimate~-df$Nc_nov)
summary(lmish)
cor.test(df$trans_Estimate,-df$Nc_nov, method = 'spearman')
cor.test(df$trans_Estimate,-df$Nc_nov, method = 'pearson')

plotISDNctrans <- ggplot(data= df, aes(-Nc_nov,trans_Estimate)) + geom_point(size=2) + labs( x = "Codon Adaptation [-Nc]", y = "Effect on ISD \n in Transmembrane Domains") 
plotISDNctrans <- plotISDNctrans + stat_smooth(method = "lm", color="red")+theme(axis.title=element_text(size=24, face ="bold"),title=element_text(size=24, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
plotISDNctrans <- plotISDNctrans+geom_text(data = df, aes(x = -62.1, y = -0.0065, label = "slope: 0.00496 \n std error: 0.00025 \n Spearman's R:0.52\n p-value:3.74e-09"),size=6.5, colour = 'black')+ theme(axis.text = element_text(size = 24)) 
plotISDNctrans


lmish<-lm(df$nontrans_Estimate~-df$Nc_nov, weights=(1/(df$nontrans_Std._Error))^2)
summary(lmish)
cor.test(df$nontrans_Estimate,-df$Nc_nov, method = 'spearman')
cor.test(df$nontrans_Estimate,-df$Nc_nov, method = 'pearson')

plotISDNcnontrans <- ggplot(data= df, aes(-Nc_nov,nontrans_Estimate)) + geom_point(size=2) + labs( x = "Codon Adaptation [-Nc]", y = "Effect on ISD \n in Nontransmembrane Pfams") 
plotISDNcnontrans <- plotISDNcnontrans + stat_smooth(method = "lm", color="red")+theme(axis.title=element_text(size=24, face ="bold"),title=element_text(size=24, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
plotISDNcnontrans <- plotISDNcnontrans +geom_text(data = df, aes(x = -62.1, y = -0.019, label = "slope: 0.00643 \n std error: 0.00049 \n Spearman's R:0.62 \n p-value:7.1e-14"),size=6.5, colour = 'black')+ theme(axis.text = element_text(size = 24)) 
plotISDNcnontrans

##################PIC################

#######read in data (ISD, clustering, etc)

namevector<-as.vector(df$SpeciesName)#vector of all the vertabrate species names

droptree<-drop.tip(tree,namevector,trim.internal=TRUE, rooted=TRUE)#dropping the tips with the vertebrate names

trimmedtree<-keep.tip(tree,namevector)# keeping only the tips with vertebrate species names


#visualize the trees
pr.species<-c("Saccharomyces_cerevisiae","Ciona_intestinalis")
tree.noPR<-drop.tip(trimmedtree,pr.species)


plotTree(tree.noPR,type="fan",fsize=0.7,lwd=1,ftype="i")
plotTree(tree.noPR,fsize = 0.4,ftype="i")



#Set up data into vectors for use in PIC
GCvector<-as.vector(df$Genomic_Percent_GC)
CAAIvector<-as.vector(df$fi_CAAI)
CAIvector<-as.vector(df$CAI)
ENCvector<-as.vector(-df$Nc_nov)
ISDvectortrans<-as.vector(df$trans_Estimate)
ISDvectornontrans<-as.vector(df$nontrans_Estimate)


######name ship

names(GCvector)<-names(CAAIvector)<-names(ENCvector)<-names(CAIvector)<- names(ISDvectortrans)<-names(ISDvectornontrans)<-namevector


######################### pic comparisons, note that the -1 forces the comparisons through the origin
piccommandGC<-pic(GCvector,trimmedtree)
piccommandCAI<- pic(CAIvector,trimmedtree)
piccommandCAAI<- pic(CAAIvector,trimmedtree)
piccommandISDtrans<-pic(ISDvectortrans, trimmedtree)
piccommandISDnontrans<-pic(ISDvectornontrans, trimmedtree)
piccommandENC<-pic(ENCvector,trimmedtree)


###############3linear models
#trans
confound_controlledISDtrans <-(lm(piccommandISDtrans~ piccommandENC-1))
summary(confound_controlledISDtrans)
cor.test(piccommandISDtrans, piccommandENC,  method = "pearson")
cor.test(piccommandISDtrans, piccommandENC,  method = "spearman")

plot<-ggplot(data= confound_controlledISDtrans, aes(piccommandENC,piccommandISDtrans)) + geom_point(size=2) + labs(x="Codon Adaptation\n[-Nc,PIC Corrected]", y="Effect on ISD [PIC Corrected] \n in Transmembrane Pfams")
plot <- plot + theme(axis.title=element_text(size=24, face ="bold"),title=element_text(size=24, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
plot<- plot+geom_text(data = df, aes(x = 0.03, y = -0.003, label = "Spearman's R:0.0033\n p-value:0.97"), colour = 'black', size = 6.5)+ theme(axis.text = element_text(size = 24)) 
plot

#nontrans
confound_controlledISDnontrans <-(lm(piccommandISDnontrans~ piccommandENC-1))
summary(confound_controlledISDnontrans)
cor.test(piccommandISDnontrans, piccommandENC,  method = "pearson")
cor.test(piccommandISDnontrans, piccommandENC,  method = "spearman")

plot1<-ggplot(data= confound_controlledISDnontrans, aes(piccommandENC,piccommandISDnontrans)) + geom_point(size=2) + labs(x="Codon Adaptation\n[-Nc,PIC Corrected]", y="Effect on ISD [PIC Corrected] \nin Nontransmembrane Pfams")
plot1 <- plot1 + theme(axis.title=element_text(size=24, face ="bold"),title=element_text(size=24, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
plot1<- plot1+geom_text(data = df, aes(x = 0.03, y = -0.021, label = "Spearman's R:0.44\n p-value:1.12e-06"), colour = 'black', size = 6.5)+ theme(axis.text = element_text(size = 24)) 
plot1


#figures
figure1 <- list(plotISDNctrans, plot)
plot_layout <- rbind(c(1,2)) # specifying a grid layout- 1 will be twice as wide as 2
grid.arrange(grobs=figure1, layout_matrix = plot_layout) 

figure2<-list(plotISDNcnontrans,plot1)
plot_layout <- rbind(c(1,2)) # specifying a grid layout- 1 will be twice as wide as 2
grid.arrange(grobs=figure2, layout_matrix = plot_layout) 

lay <- rbind(c(1,2,3),
             c(1,4,5))

grid.arrange(ISD3,plotISDNctrans, plot,plotISDNcnontrans,plot1, layout_matrix = lay)

