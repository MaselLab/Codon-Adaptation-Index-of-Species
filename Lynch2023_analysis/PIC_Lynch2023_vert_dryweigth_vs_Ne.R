#For my analysis I followed https://www.zoology.ubc.ca/~schluter/R/Phylogenetic.html
library(xlsx)
library(ggplot2)
library(ape)
library(rstudioapi)

setwd(dirname(getSourceEditorContext()$path))

longlist <- read.xlsx(file = 'Data/Lynch2023_embr202357561-sup-0002-datasetev1.xlsx', 1)

beg_vert = which(longlist$Dataset.EV1 == 'Vertebrates:')
end_vert = which(longlist$Dataset.EV1 == 'Vascular plants:')

vert_list <- longlist[(beg_vert+1):(end_vert-1), ]
colnames(vert_list) <- longlist[14, ]
#Clean spaces between columns and rows
vert_list <- vert_list[ ,colSums(is.na(vert_list)) < nrow(vert_list)]
vert_list <- vert_list[rowSums(is.na(vert_list)) < ncol(vert_list), ]

vert_ana <- data.frame(Species = vert_list[,1],
                         logMutation_rate = log10(as.numeric(vert_list[,2])), 
                         logDry_weight = log10(as.numeric(vert_list[,12])),
                         logNe = log10(as.numeric(vert_list$Ne)))

#Some species have multiple values for mu
#I'm only taken the first value reported for any mu. Take mean instead?
vert_ana <- na.omit(vert_ana)

ggplot(vert_ana, aes(logDry_weight, logNe))+ 
  geom_point()

model <- lm(logNe ~ logDry_weight, data = vert_ana)
model
summary(model)

#Cleaning some of the names in the database to match the database of timetree
vert_ana$Species[which(vert_ana$Species == 'Sus scrofa / suis')] = 'Sus scrofa'
vert_ana$Species[which(vert_ana$Species == 'Saimiri boliviensis boliviensis')] = 'Saimiri boliviensis'

write(vert_ana$Species, 'Data/vert_species_indata.txt')
#I create the tree by uploading manually Data/vert_species_indata.txt into http://timetree.org/
#After this vert_species_indata.nwk is created in the Download folder
#Manually I move this file to the Data folder
vert_tree <- read.tree("Data/vert_species_indata.nwk") 
plot(vert_tree, cex=0.7)

#Usage of name replacements to set the same Species name in the data as in the tree
#Note that when TimeTree replaces spaces by _
#Last command does not deal with Species with two spaces
#Because they are not too many I deal with this manually
vert_ana$Species <- sub(" ", "_", vert_ana$Species)
vert_ana$Species[which(vert_ana$Species == 'Cervus_elaphus yarkandensis')] = 'Cervus_elaphus_yarkandensis'
vert_ana$Species[which(vert_ana$Species == 'Ceratotherium_simum simum')] = 'Ceratotherium_simum_simum'
vert_ana$Species[which(vert_ana$Species == 'Saimiri_boliviensis boliviensis')] = 'Saimiri_boliviensis_boliviensis'
#Neovison_vison is recognized in Time Tree as Neogale_vison
vert_ana$Species[which(vert_ana$Species == 'Neovison_vison')] = 'Neogale_vison'

rownames(vert_ana) <- vert_ana$Species

#Check that all the species listed on the data frame are in the tree and all 
#species listed on the tree are on the data frame
vert_ana$Species[which(is.na(match(rownames(vert_ana), vert_tree$tip.label)))]
vert_tree$tip[which(is.na(match(vert_tree$tip.label, rownames(vert_ana))))]

vert_ana <- vert_ana[match(vert_tree$tip.label, rownames(vert_ana)),]

#Calculate pic corrected variables
pic_vert_ana <- data.frame(PIC_logDry_weight <- pic(vert_ana$logDry_weight, vert_tree),
                           PIC_logNe <- pic(vert_ana$logNe, vert_tree))

colnames(pic_vert_ana) <- c('PIC_logDry_weight', 'PIC_logNe')

ggplot(pic_vert_ana, aes(PIC_logDry_weight, PIC_logNe))+ 
  geom_point()

model <- lm(PIC_logNe ~ PIC_logDry_weight, data = pic_vert_ana)
model
summary(model)

#export pic data

write.csv(pic_vert_ana, 'Data/pic_vert.txt')

