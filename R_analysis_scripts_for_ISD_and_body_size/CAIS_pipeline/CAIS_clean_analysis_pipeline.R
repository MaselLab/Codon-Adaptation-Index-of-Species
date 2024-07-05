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

#merge
df<-merge(dfvert, dfGC, by.x="SpeciesName", by.y="NewickName")
str(df)

######################RAW############

hist(df$Estimate)

lmish<-lm(df$Estimate~df$fi_CAAI)
summary(lmish)

weightCAIS<-lm(df$Estimate~df$fi_CAAI,weights=(1/(df$Std._Error))^(2))
summary(weightCAIS)

cor.test(df$Estimate,df$fi_CAAI, method = 'pearson')
cor.test(df$Estimate,df$fi_CAAI, method = 'spearman')

lmGC<-lm(df$Genomic_Percent_GC~df$fi_CAAI)
summary(lmGC)
cor.test(df$Genomic_Percent_GC,df$fi_CAAI, method='spearman')

plotISDCAI <- ggplot(data= df, aes(fi_CAAI,Estimate)) + geom_point(size=2) + labs( x = "CAIS", y = "Effect on ISD \n in Young Domains") 
plotISDCAI <- plotISDCAI + stat_smooth(method = "lm", color="red")+theme(axis.title=element_text(size=28, face ="bold"),title=element_text(size=24, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
plotISDCAI <- plotISDCAI+geom_text(data = df, aes(x = 1.061, y = 0.02, label = "Spearman's R:0.57  \n Pearson's R:0.70\n p-value:2e-11"),size=7, colour = 'black')+ theme(axis.text = element_text(size = 24)) 
plotISDCAI


plotISDGC <-ggplot(data= df, aes(Genomic_Percent_GC,Estimate, color=fi_CAAI)) + geom_point(size=2) + labs(title= "% GC", x = "Genomic Percent GC", y = "Effect on ISD ") 
plotISDGC <- plotISDGC + stat_smooth(method = "lm", color="red")+ theme(axis.text = element_text(size = 14)) 
plotISDGC <-plotISDGC +geom_text(data = df, aes(x = 0.375, y = 0.0075, label = "Spearman's R: 0.073 \n Pearson's R:-0.067\n p-value:0.07"), colour = 'black')
plotISDGC

####GC _ CAI


plotGCCAI <- ggplot(data= df, aes(fi_CAAI,Genomic_Percent_GC)) + geom_point(size=2)+ labs( y = "Genomic Percent GC", x = "Codon Adaptation Index of Species (CAIs)")+theme(axis.title=element_text(size=11, face ="bold", size = 14),title=element_text(size=14, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
plotGCCAI <-plotGCCAI + geom_text(data = df, aes(x = 1.065, y = 0.44, label = "Spearman's R: 0.073 \n Pearson's R:-0.067\n p-value:0.48"), colour = 'black')+ theme(axis.text = element_text(size = 14)) 
plotGCCAI

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

######name ship

names(GCvector)<-names(CAAIvector)<-names(CAIvector)<- names(ISDvector)<-namevector


######################### pic comparisons, note that the -1 forces the comparisons through the origin
piccommandGC<-pic(GCvector,trimmedtree)
piccommandCAI<- pic(CAIvector,trimmedtree)
piccommandCAAI<- pic(CAAIvector,trimmedtree)
piccommandISD<-pic(ISDvector, trimmedtree)


###############3linear models
confound_controlledISD <-(lm(piccommandISD~ piccommandCAAI-1))
summary(confound_controlledISD)
cor.test(piccommandISD, piccommandCAAI,  method = "pearson")
cor.test(piccommandISD, piccommandCAAI,  method = "spearman")

confound_controlledGC <-(lm(piccommandCAAI~ piccommandGC-1))
summary(confound_controlledGC)
cor.test(piccommandGC, piccommandCAAI,  method = "pearson")

hist(piccommandCAAI)
hist(piccommandISD)



plot2 <-ggplot(data= confound_controlledGC, aes(piccommandCAAI,piccommandGC)) + geom_point(size=2) + labs( x = "CAIS [PIC Corrected]", y = "%GC [PIC Corrected]") 
plot2 <- plot2+theme(axis.title=element_text(size=28, face ="bold"),title=element_text(size=28, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
plot2 <- plot2 +geom_text(data = df, aes(x = 0.0022, y = -0.002, label = "Spearman's R: 0.23 \n Pearson's R:0.43\n p-value:0.012"), colour = 'black', size = 7)+ theme(axis.text = element_text(size = 24)) 
plot2

plot<-ggplot(data= confound_controlledISD, aes(piccommandCAAI,piccommandISD)) + geom_point(size=2) + labs(x="CAIS [PIC Corrected]", y="Effect on ISD [PIC Corrected]\n in Young Domains")+geom_smooth(method = "lm", color="red")
plot <- plot + theme(axis.title=element_text(size=28, face ="bold"),title=element_text(size=14, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
plot<- plot+geom_text(data = df, aes(x = -0.0014, y = 0.0011, label = "Spearman's R: 0.36 \n Pearson's R:0.47\n p-value:2e-7"), colour = 'black', size = 7)+ theme(axis.text = element_text(size = 24)) 
plot

plot1<-ggplot(data= confound_controlledISD, aes(piccommandGC,piccommandISD, color=piccommandCAAI)) + geom_point(size=2) + labs(title="Phylogenetically Corrected Genomic Percent GC \n vs \n Phylogenetically Corrected Species Effect on ISD", x="Phylogenetically Corrected Genomic Percent GC", y="Phylogenetically Corrected Species Effect on ISD")+geom_smooth(method = "lm", color="red")+ theme(axis.text = element_text(size = 14)) 
plot1<- plot1+geom_text(data = df, aes(x = -0.001, y = 0.0017, label = "Spearman's R:0.37 \n Pearson's R:0.41\n p-value:4.8e-06"), colour = 'black')+theme(axis.title=element_text(size=12, face ="bold"),title=element_text(size=14, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
plot1

plot3 <-ggplot(data= confound_controlledGC, aes(piccommandCAI,piccommandGC)) + geom_point(size=2) + labs( x = "CAI [PIC Corrected]", y = "%GC [PIC Corrected]") 
plot3 <- plot3+theme(axis.title=element_text(size=28, face ="bold"),title=element_text(size=28, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
plot3 <- plot3 +geom_text(data = df, aes(x = 0.006, y = 0.002, label = "Spearman's R: -0.6 \n Pearson's R:-0.55\n p-value:1.1e-10"), colour = 'black', size = 7)+ theme(axis.text = element_text(size = 24)) 
plot3

#figures
figure1 <- list(plotISDCAI, plot)
plot_layout <- rbind(c(1,2)) # specifying a grid layout- 1 will be twice as wide as 2
grid.arrange(grobs=figure1, layout_matrix = plot_layout) 

figure2<-list(plotGCCAI,plot2)
plot_layout <- rbind(c(1,2)) # specifying a grid layout- 1 will be twice as wide as 2
grid.arrange(grobs=figure2, layout_matrix = plot_layout) 

figure4<-list(massplot3,plot2,massplot2,plot3)
#plot_layout <- rbind(c(1,2)) # specifying a grid layout- 1 will be twice as wide as 2
#grid.arrange(grobs=figure4, layout_matrix = plot_layout)

grid.arrange(CAI_CAIMAX,plot3,plot2,CAI_CAAI2,massplot3,massplot2, nrow=2)



