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

dfmass<-read.csv(file = "body_mass.txt", header = T)
str(dfmass)

dfno_weights<-read.csv(file="CAIS_table.txt", header=T)
str(dfno_weights)

dfspecies<-read.csv(file="SpeciesList_11_7_2019.txt",header=T)
str(dfspecies)

dfoldCAI<-read.csv(file = "GC_CAI_data.txt", header = T)
str(dfoldCAI)

dfintergenic<-read.csv(file = "Intergenic_CAIS.txt", header= T)
str(dfintergenic)

dfENC<-read.csv(file = "ENC_6_8_2020.txt",header = T)
str(dfENC)

dfgeneGC<-read.csv(file = "Species_Gene_GC_data.txt",header = T)
str(dfgeneGC)
dfmass<-merge(dfmass, dfspecies, by.x="NewickName", by.y="NewickName")
str(dfmass)

dfmass<-merge(dfmass, dfENC, by.x="SpeciesUID", by.y="SpeciesUID")
str(dfmass)

dfmass<-merge(dfmass, dfno_weights, by.x="SpeciesUID", by.y="SpeciesUID")
str(dfmass)

dfmass<-merge(dfmass, dfintergenic, by.x="SpeciesUID", by.y="SpeciesUID")
str(dfmass)

dfmass<-merge(dfmass, dfgeneGC, by.x="SpeciesUID", by.y="SpeciesUID")
str(dfmass)

#read in newick tree
tree<-read.tree(file = "PhylogeneticTree_AllSpecies.nwk.txt")

lmi<-lm(dfmass$X5.1_AdultBodyMass_g~log10(dfmass$complete_CAIS.x))
summary(lmi)
cor.test(dfmass$complete_CAIS.x,log10(dfmass$X5.1_AdultBodyMass_g), method = 'pearson')
cor.test(dfmass$complete_CAIS.x,log10(dfmass$X5.1_AdultBodyMass_g), method = 'spearman')

cor.test(dfmass$complete_CAIS.x,dfmass$X5.1_AdultBodyMass_g, method = 'pearson')
cor.test(dfmass$complete_CAIS.x,dfmass$X5.1_AdultBodyMass_g, method = 'spearman')

massplot <-ggplot(data= dfmass, aes(complete_CAIS.x,log10(X5.1_AdultBodyMass_g))) + geom_point(size=2) + labs( x = "Codon Adaptation Index of Species", y = "Adult Body Mass (g) (log 10)") + stat_smooth(method = "lm", color="red") +stat_smooth(method = "loess", color="green")
massplot <- massplot+theme(axis.title=element_text(size=14, face ="bold"),title=element_text(size=14, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
massplot

massplot1 <-ggplot(data= dfmass, aes(complete_CAIS.x,X5.1_AdultBodyMass_g)) + geom_point(size=2) + labs( x = "Codon Adaptation Index of Species", y = "Adult Body Mass (g)") + stat_smooth(method = "lm", color="red") +stat_smooth(method = "loess", color="green")
massplot1 <- massplot1+theme(axis.title=element_text(size=14, face ="bold"),title=element_text(size=14, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
massplot1

###################LocalGC########################### something funky this way grows
cor.test(dfmass$complete_CAIS.x,dfmass$gene_GC_variance, method = 'pearson')
cor.test(dfmass$complete_CAIS.x,dfmass$gene_GC_variance, method = 'spearman')

plot(x=dfmass$complete_CAIS.x,y=dfmass$gene_GC_variance)

cor.test(dfmass$CAIS_intergenic,dfmass$gene_GC_variance, method = 'pearson')
cor.test(dfmass$CAIS_intergenic,dfmass$gene_GC_variance, method = 'spearman')

plot(x=dfmass$CAIS_intergenic,y=dfmass$gene_GC_variance)

##################PIC###############

#######read in data (ISD, clustering, etc)

namevector<-as.vector(dfmass$NewickName.x)#vector of all the vertabrate species names

droptree<-drop.tip(tree,namevector,trim.internal=TRUE, rooted=TRUE)#dropping the tips with the vertebrate names

trimmedtree<-keep.tip(tree,namevector)# keeping only the tips with vertebrate species names


#visualize the trees
pr.species<-c("Saccharomyces_cerevisiae","Ciona_intestinalis")
tree.noPR<-drop.tip(trimmedtree,pr.species)


plotTree(tree.noPR,type="fan",fsize=0.7,lwd=1,ftype="i")
plotTree(tree.noPR,fsize = 0.4,ftype="i")

#Set up data into vectors for use in PIC
massvector<-as.vector(log10(dfmass$X5.1_AdultBodyMass_g))
CAISvector<-as.vector(dfmass$complete_CAIS.x)
#CAISvector<-as.vector(dfmass$Percent_GC_intergenic)
CAIvector<-as.vector(dfmass$CAI.y)
ENCvector<-as.vector(dfmass$Nc_nov)
equalCAISvector<-as.vector(dfmass$CAIS_no_AA)
GCvector<-as.vector(dfmass$Genomic_Percent_GC.x)

######name ship

names(CAISvector)<-names(massvector)<-names(CAIvector)<-names(ENCvector)<-names(equalCAISvector)<-names(GCvector)


######################### pic comparisons, note that the -1 forces the comparisons through the origin
piccommandmass<-pic(massvector,trimmedtree)
piccommandCAIS<- pic(CAISvector,trimmedtree)
piccommandCAI<- pic(CAIvector,trimmedtree)
piccommandENC<- pic(ENCvector,trimmedtree)
piccommandequalCAIS<- pic(equalCAISvector,trimmedtree)
piccommandGC<- pic(GCvector,trimmedtree)

lmiPIC<-lm(piccommandmass~ piccommandCAIS-1)
summary(lmiPIC)
cor.test(piccommandmass, piccommandCAIS, method = 'pearson')
cor.test(piccommandmass, piccommandCAIS, method = 'spearman')

massplot2 <-ggplot(data=lmiPIC, aes(piccommandCAIS,piccommandmass)) + geom_point(size=2) + labs( x = "CAIS [PIC Corrected]", y = "Adult Body Mass\n[PIC Corrected]") + stat_smooth(method = "lm", color="red")
massplot2 <- massplot2+theme(axis.title=element_text(size=28, face ="bold"),title=element_text(size=18, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
massplot2<- massplot2 + geom_text(data=lmiPIC, aes(x = -0.0005, y = -0.8, label = "Pearson's R:-0.33\n p-value:0.0065"),size=7, colour = 'black')+ theme(axis.text = element_text(size = 24)) 
massplot2

massplot2a <-ggplot(data=lmiPIC, aes(piccommandCAIS,piccommandmass)) + geom_point(size=2) + labs( x = "CAIS [PIC Corrected] \n Calculated from Intergenic % GC", y = "Adult Body Mass\n[PIC Corrected]") + stat_smooth(method = "lm", color="red")
massplot2a <- massplot2a+theme(axis.title=element_text(size=28, face ="bold"),title=element_text(size=18, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
massplot2a<- massplot2a + geom_text(data=lmiPIC, aes(x = -0.005, y = -0.8, label = "Pearson's R:-0.46\n p-value:0.00014"),size=7, colour = 'black')+ theme(axis.text = element_text(size = 24)) 
massplot2a


lmiPIC1<-lm(piccommandmass~ piccommandCAI-1)
summary(lmiPIC1)
cor.test(piccommandmass, piccommandCAI, method = 'pearson')
cor.test(piccommandmass, piccommandCAI, method = 'spearman')

massplot3 <-ggplot(data=lmiPIC1, aes(piccommandCAI,piccommandmass)) + geom_point(size=2) + labs( x = "CAI [PIC Corrected]", y = "Adult Body Mass\n[PIC Corrected]") + stat_smooth(method = "lm", color="red")
massplot3 <- massplot3+theme(axis.title=element_text(size=28, face ="bold"),title=element_text(size=18, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
massplot3 <- massplot3+geom_text(data=lmiPIC, aes(x = 0.002, y = -0.8, label = "Pearson's R:0.297\n p-value:0.015"),size=7, colour = 'black')+ theme(axis.text = element_text(size = 24)) 
massplot3

lmiPIC<-lm(piccommandmass~ piccommandENC-1)
summary(lmiPIC)
cor.test(piccommandmass, piccommandENC, method = 'pearson')
cor.test(piccommandmass, piccommandENC, method = 'spearman')

massplot4 <-ggplot(data=lmiPIC, aes(piccommandENC,piccommandmass)) + geom_point(size=2) + labs( x = "ENC [PIC Corrected]", y = "Adult Body Mass\n[PIC Corrected]") + stat_smooth(method = "lm", color="red")
massplot4 <- massplot4+theme(axis.title=element_text(size=28, face ="bold"),title=element_text(size=18, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
massplot4<- massplot4+geom_text(data=lmiPIC, aes(x = 0.002, y = -0.8, label = "Pearson's R:-0.00825\n p-value:0.948"),size=7, colour = 'black')+ theme(axis.text = element_text(size = 24)) 
massplot4

lmiPIC<-lm(piccommandmass~ piccommandequalCAIS-1)
summary(lmiPIC)
cor.test(piccommandmass, piccommandequalCAIS, method = 'pearson')
cor.test(piccommandmass, piccommandequalCAIS, method = 'spearman')

massplot5 <-ggplot(data=lmiPIC, aes(piccommandequalCAIS,piccommandmass)) + geom_point(size=2) + labs( x = "CAIS without AA Correction \n [PIC Corrected]", y = "Adult Body Mass\n[PIC Corrected]") + stat_smooth(method = "lm", color="red")
massplot5 <- massplot5+theme(axis.title=element_text(size=28, face ="bold"),title=element_text(size=18, face ="bold"),panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey90"))
massplot5<- massplot5+geom_text(data=lmiPIC, aes(x = 0.002, y = -0.8, label = "Pearson's R:-0.348\n p-value:0.0044"),size=7, colour = 'black')+ theme(axis.text = element_text(size = 24)) 
massplot5

cor.test(piccommandmass, piccommandGC, method = 'pearson')
cor.test(piccommandmass, piccommandGC, method = 'spearman')
cor.test(piccommandequalCAIS, piccommandGC, method = 'pearson')
cor.test(piccommandequalCAIS, piccommandGC, method = 'spearman')
cor.test(piccommandCAIS, piccommandGC, method = 'pearson')
cor.test(piccommandCAIS, piccommandGC, method = 'spearman')
cor.test(piccommandENC, piccommandGC, method = 'pearson')
cor.test(piccommandENC, piccommandGC, method = 'spearman')

figure<-list(massplot3,massplot2,massplot4)
plot_layout <- rbind(c(1,3)) # specifying a grid layout- 1 will be twice as wide as 2
grid.arrange(grobs=figure, layout_matrix = plot_layout) 

grid.arrange(massplot2,massplot3,massplot4, nrow=1)

