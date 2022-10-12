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
library(stats)

#read in newick tree
tree<-read.tree(file = "PhylogeneticTree_AllSpecies.nwk.txt")


####read in pertinent data
##########
#dfvert<-read.csv(file = "clean_ISD_effects.txt", header = T)
#str(dfvert)

#dfENC<-read.csv(file = "ENC_6_8_2020.txt",header = T)
#str(dfENC)

#dfvert<-read.csv(file = "clean_recent_ISD_effects.txt", header = T)
#str(dfvert)

dfvert<-read.csv(file = "clean_ancient_ISD_effects.txt", header = T)
str(dfvert)

#clean out unneccesary bits
dfvert<-subset(dfvert, select = -c(df))
str(dfvert)


dfGC<-read.csv(file = "CAAI_CAI_GC_1_29_2020.txt", header = T)
str(dfGC)

dfGC<-subset(dfGC, select = -c(Row,Speciesname))
str(dfGC)

dfGC<-merge(dfGC, dfENC, by.x="SpeciesUID", by.y="SpeciesUID")
str(dfGC)

#merge
df<-merge(dfvert, dfGC, by.x="SpeciesName", by.y="NewickName")
str(df)


######################RAW#############

hist(df$Estimate)

lmish<-lm(df$Estimate~-df$Nc_nov)
summary(lmish)

weightENC<-lm(df$Estimate~-df$Nc_nov,weights=(1/(df$Std._Error))^(2))
summary(weightENC)

cor.test(df$Estimate,-df$Nc_nov, method = 'pearson')
cor.test(df$Estimate,-df$Nc_nov, method = 'spearman')

lmGC<-lm(df$Genomic_Percent_GC~df$Nc_nov)
summary(lmGC)
cor.test(df$Genomic_Percent_GC,df$Nc_nov, method='pearson')

plotISDNc <- ggplot(data= df, aes(-Nc_nov,Estimate)) + geom_point(size=2) + labs( x = "Codon Adaptation [-Nc]", y = "Effect on ISD \n in Old Domains") 
plotISDNc <- plotISDNc + stat_smooth(method = "lm", color="red")+theme(axis.title=element_text(size=24, face ="bold"),axis.text=element_text(size=24, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
plotISDNc <- plotISDNc+geom_text(data = df, aes(x = -62.85, y = 0.015, label = "slope: 0.0073\n std error:0.0003089 \n Spearman's R:0.64\n Pearson's R:0.69\n p-value:1.24e-14"),size=6.5, colour = 'black')+ theme(axis.text = element_text(size = 24)) 
plotISDNc

####GC _ Nc


plotGCnc <- ggplot(data= df, aes(Nc_nov,Genomic_Percent_GC)) + geom_point(size=2)+ labs( y = "%GC", x = "Nc")+theme(axis.title=element_text(size=24, face ="bold"),axis.text = element_text(size=24, face= "bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
plotGCnc <-plotGCnc + geom_text(data = df, aes(x = 62, y = 0.38, label = "Spearman's R: 0.18 \n Pearson's R:0.29\n p-value:0.05229"), size=6.5, colour = 'black') 
plotGCnc

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
ISDvector<-as.vector(df$Estimate)
ENCvector<-as.vector(-df$Nc_nov)
######name ship

names(GCvector)<-names(CAAIvector)<-names(CAIvector)<- names(ISDvector)<-names(ENCvector)<-namevector


######################### pic comparisons, note that the -1 forces the comparisons through the origin
piccommandGC<-pic(GCvector,trimmedtree)
piccommandCAI<- pic(CAIvector,trimmedtree)
piccommandCAAI<- pic(CAAIvector,trimmedtree)
piccommandISD<-pic(ISDvector, trimmedtree)
piccommandENC<-pic(ENCvector,trimmedtree)

###############3linear models
confound_controlledISD <-(lm(piccommandISD~ piccommandCAAI-1))
summary(confound_controlledISD)
cor.test(piccommandISD, piccommandCAAI,  method = "pearson")
cor.test(piccommandISD, piccommandCAAI,  method = "spearman")

confound_controlledGC <-(lm(piccommandENC~ piccommandGC-1))
summary(confound_controlledGC)
cor.test(piccommandGC, piccommandENC,  method = "pearson")

confound_controlledENC <-(lm(piccommandISD~ piccommandENC-1))
summary(confound_controlledENC)
cor.test(piccommandISD, piccommandENC,  method = "pearson")
cor.test(piccommandISD, piccommandENC,  method = "spearman")


plot2 <-ggplot(data= confound_controlledGC, aes(piccommandENC,piccommandGC)) + geom_point(size=2) + labs( x = "Nc [PIC Corrected]", y = "%GC [PIC Corrected]") 
plot2 <- plot2+theme(axis.title=element_text(size=24, face ="bold"),title=element_text(size=24, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
plot2 <- plot2 +geom_text(data = df, aes(x = -0.04, y = -0.002, label = "Spearman's R: 0.099\n Pearson's R:0.12\n p-value:0.29"), colour = 'black', size = 6.5)+ theme(axis.text = element_text(size = 24)) 
plot2

plot<-ggplot(data= confound_controlledENC, aes(piccommandENC,piccommandISD)) + geom_point(size=2) + labs(x="Codon Adaptation [-Nc,PIC Corrected]", y="Effect on ISD [PIC Corrected] \n in Old Domains")+geom_smooth(method = "lm", color="red")
plot <- plot + theme(axis.title=element_text(size=24, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
plot<- plot+geom_text(data = df, aes(x = -0.031, y = 0.0008, label = "Spearman's R: 0.43 \n Pearson's R:0.44\n p-value:1.936e-06"), colour = 'black', size = 6.5)+ theme(axis.text = element_text(size = 24)) 
plot


#figures
figure1 <- list(plotISDNc, plot)
plot_layout <- rbind(c(1,2)) # specifying a grid layout- 1 will be twice as wide as 2
grid.arrange(grobs=figure1, layout_matrix = plot_layout) 

figure2 <- list(plotGCnc, plot2)
plot_layout <- rbind(c(1,2)) # specifying a grid layout- 1 will be twice as wide as 2
grid.arrange(grobs=figure2, layout_matrix = plot_layout)


#grid.arrange(CAI_CAIMAX,plot3,plot2,CAI_CAAI2,massplot3,massplot2, nrow=2)



