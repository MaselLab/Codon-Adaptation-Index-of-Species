library(MASS)
library(nlme)
library(lme4)
require(lmerTest)
library(optimx)
library(ggplot2)
library(MuMIn)
library(glmm)

####read in pertinent data
dfNCBI<-read.csv(file = "ISD_clust_dataNCBI_vertebrates.txt", header = T)
dfNCBI$ISD_WithCys <-as.character((dfNCBI$ISD_WithCys))
dfNCBI$ISD_WithCys <-as.numeric((dfNCBI$ISD_WithCys))
dfNCBI$Clustering <-as.character((dfNCBI$Clustering))
dfNCBI$Clustering <-as.numeric((dfNCBI$Clustering))
dfNCBI$TmhmmTopology <-as.character((dfNCBI$TmhmmTopology))
dfNCBI$TmhmmTopology <-as.numeric((dfNCBI$TmhmmTopology))
dfNCBI$SpeciesUID <-as.character((dfNCBI$SpeciesUID))
dfNCBI$SpeciesUID <-as.integer((dfNCBI$SpeciesUID))
dfNCBI$UID <-as.character((dfNCBI$UID))
dfNCBI$UID <-as.integer((dfNCBI$UID))
str(dfNCBI)

dfensembl<-read.csv(file = "ISD_clust_dataensembl_vertebrates.txt", header = T)
dfensembl$ISD_WithCys <-as.character((dfensembl$ISD_WithCys))
dfensembl$ISD_WithCys <-as.numeric((dfensembl$ISD_WithCys))
dfensembl$Clustering <-as.character((dfensembl$Clustering))
dfensembl$Clustering <-as.numeric((dfensembl$Clustering))
dfensembl$TmhmmTopology <-as.character((dfensembl$TmhmmTopology))
dfensembl$TmhmmTopology <-as.numeric((dfensembl$TmhmmTopology))
str(dfensembl)

total <- rbind(dfNCBI, dfensembl)

dfSpecies<-read.csv(file = "SpeciesList_11_7_2019.txt", header = T)
dfSpecies<-subset(dfSpecies, select = c(SpeciesUID,NewickName))

dfages<-read.csv(file = "~/LMM_project/Pfamages.txt", header = T)

total<-merge(total, dfages, by.x="PfamUID", by.y="PfamUID")
total<-merge(total, dfSpecies, by.x="SpeciesUID", by.y="SpeciesUID")


###########

cleanvertebratelist<- c('Anas_platyrhynchos' , 'Anolis_carolinensis' , 'Anser_cygnoides' , 'Apaloderma_vittatum' , 'Balaenoptera_acutorostrata' , 'Bison_bison' , 'Boleophthalmus_pectinirostris' , 'Bos_indicus' , 'Bos_mutus' , 'Bubalus_bubalis' , 'Buceros_rhinoceros' , 'Calidris_pugnax' , 'Callithrix_jacchus' , 'Callorhinchus_milii' , 'Calypte_anna' , 'Camelus_bactrianus' , 'Camelus_dromedarius' , 'Camelus_ferus' , 'Cariama_cristata' , 'Carlito_syrichta' , 'Castor_canadensis' , 'Cavia_porcellus' , 'Cebus_capucinus' , 'Ceratotherium_simum' , 'Cercocebus_atys' , 'Charadrius_vociferus' , 'Chelonia_mydas' , 'Chlorocebus_sabaeus' , 'Chrysemys_picta' , 'Chrysochloris_asiatica' , 'Clupea_harengus' , 'Colobus_angolensis' , 'Columba_livia' , 'Condylura_cristata' , 'Corvus_brachyrhynchos' , 'Coturnix_japonica' , 'Cricetulus_griseus' , 'Crocodylus_porosus' , 'Cuculus_canorus' , 'Cyprinodon_variegatus' , 'Cyprinus_carpio' , 'Danio_rerio' , 'Dipodomys_ordii' , 'Elephantulus_edwardii' , 'Eptesicus_fuscus' , 'Equus_caballus' , 'Equus_przewalskii' , 'Erinaceus_europaeus' , 'Esox_lucius' , 'Eurypyga_helias' , 'Falco_peregrinus' , 'Felis_catus' , 'Fukomys_damarensis' , 'Fulmarus_glacialis' , 'Fundulus_heteroclitus' , 'Gekko_japonicus' , 'Geospiza_fortis' , 'Gorilla_gorilla' , 'Hzliaeetus_albicilla' , 'Haliaeetus_leucocephalus' , 'Heterocephalus_glaber' , 'Hipposideros_armiger' , 'Homo_sapiens' , 'Ictalurus_punctatus' , 'Jaculus_jaculus' , 'Labrus_bergylta' , 'Latimeria_chalumnae' , 'Lepidothrix_coronata' , 'Lepisosteus_oculatus' , 'Lipotes_vexillifer' , 'Loxodonta_africana' , 'Macaca_mulatta' , 'Macaca_nemestrina' , 'Mandrillus_leucophaeus' , 'Marmota_marmota' , 'Meriones_unguiculatus' , 'Mesocricetus_auratus' , 'Microcebus_murinus' , 'Microtus_ochrogaster' , 'Miniopterus_natalensis' , 'Mus_caroli' , 'Mus_musculus' , 'Mus_pahari' , 'Mustela_putorius' , 'Myotis_brandtii' , 'Myotis_davidii' , 'Myotis_lucifugus' , 'Nannospalax_galili' , 'Nomascus_leucogenys' , 'Octodon_degus' , 'Odobenus_rosmarus' , 'Odocoileus_virginianus' , 'Orcinus_orca' , 'Ornithorhynchus_anatinus' , 'Oryctolagus_cuniculus' , 'Otolemur_garnettii' , 'Pan_paniscus' , 'Pan_troglodytes' , 'Panthera_pardus' , 'Pantholops_hodgsonii' , 'Pelodiscus_sinensis' , 'Peromyscus_maniculatus' , 'Physeter_catodon' , 'Pongo_abelii' , 'Propithecus_coquereli' , 'Pteropus_alecto' , 'Rattus_norvegicus' , 'Rhinopithecus_bieti' , 'Rousettus_aegyptiacus' , 'Saimiri_boliviensis' , 'Sorex_araneus' , 'Sus_scrofa' , 'Taeniopygia_guttata' , 'Takifugu_rubripes' , 'Trichechus_manatus' , 'Ursus_maritimus' , 'Xenopus_tropicalis' , 'Xiphophorus_maculatus')
cleananimallist<-c('Anas_platyrhynchos' , 'Anolis_carolinensis' , 'Anser_cygnoides' , 'Apaloderma_vittatum' , 'Balaenoptera_acutorostrata' , 'Bison_bison' , 'Boleophthalmus_pectinirostris' , 'Bos_indicus' , 'Bos_mutus' , 'Bubalus_bubalis' , 'Buceros_rhinoceros' , 'Calidris_pugnax' , 'Callithrix_jacchus' , 'Callorhinchus_milii' , 'Calypte_anna' , 'Camelus_bactrianus' , 'Camelus_dromedarius' , 'Camelus_ferus' , 'Cariama_cristata' , 'Carlito_syrichta' , 'Castor_canadensis' , 'Cavia_porcellus' , 'Cebus_capucinus' , 'Ceratotherium_simum' , 'Cercocebus_atys' , 'Charadrius_vociferus' , 'Chelonia_mydas' , 'Chlorocebus_sabaeus' , 'Chrysemys_picta' , 'Chrysochloris_asiatica' , 'Clupea_harengus' , 'Colobus_angolensis' , 'Columba_livia' , 'Condylura_cristata' , 'Corvus_brachyrhynchos' , 'Coturnix_japonica' , 'Cricetulus_griseus' , 'Crocodylus_porosus' , 'Cuculus_canorus' , 'Cyprinodon_variegatus' , 'Cyprinus_carpio' , 'Danio_rerio' , 'Dipodomys_ordii' , 'Elephantulus_edwardii' , 'Eptesicus_fuscus' , 'Equus_caballus' , 'Equus_przewalskii' , 'Erinaceus_europaeus' , 'Esox_lucius' , 'Eurypyga_helias' , 'Falco_peregrinus' , 'Felis_catus' , 'Fukomys_damarensis' , 'Fulmarus_glacialis' , 'Fundulus_heteroclitus' , 'Gekko_japonicus' , 'Geospiza_fortis' , 'Gorilla_gorilla' , 'Hzliaeetus_albicilla' , 'Haliaeetus_leucocephalus' , 'Heterocephalus_glaber' , 'Hipposideros_armiger' , 'Homo_sapiens' , 'Ictalurus_punctatus' , 'Jaculus_jaculus' , 'Labrus_bergylta' , 'Latimeria_chalumnae' , 'Lepidothrix_coronata' , 'Lepisosteus_oculatus' , 'Lipotes_vexillifer' , 'Loxodonta_africana' , 'Macaca_mulatta' , 'Macaca_nemestrina' , 'Mandrillus_leucophaeus' , 'Marmota_marmota' , 'Meriones_unguiculatus' , 'Mesocricetus_auratus' , 'Microcebus_murinus' , 'Microtus_ochrogaster' , 'Miniopterus_natalensis' , 'Mus_caroli' , 'Mus_musculus' , 'Mus_pahari' , 'Mustela_putorius' , 'Myotis_brandtii' , 'Myotis_davidii' , 'Myotis_lucifugus' , 'Nannospalax_galili' , 'Nomascus_leucogenys' , 'Octodon_degus' , 'Odobenus_rosmarus' , 'Odocoileus_virginianus' , 'Orcinus_orca' , 'Ornithorhynchus_anatinus' , 'Oryctolagus_cuniculus' , 'Otolemur_garnettii' , 'Pan_paniscus' , 'Pan_troglodytes' , 'Panthera_pardus' , 'Pantholops_hodgsonii' , 'Pelodiscus_sinensis' , 'Peromyscus_maniculatus' , 'Physeter_catodon' , 'Pongo_abelii' , 'Propithecus_coquereli' , 'Pteropus_alecto' , 'Rattus_norvegicus' , 'Rhinopithecus_bieti' , 'Rousettus_aegyptiacus' , 'Saimiri_boliviensis' , 'Sorex_araneus' , 'Sus_scrofa' , 'Taeniopygia_guttata' , 'Takifugu_rubripes' , 'Trichechus_manatus' , 'Ursus_maritimus' , 'Xenopus_tropicalis' , 'Xiphophorus_maculatus,Acyrthosiphon_pisum','Aedes_aegypti','Amphimedon_queenslandica','Anopheles_darlingi','Anopheles_gambiae','Apis_mellifera','Atta_cephalotes','Belgica_antarctica','Bombyx_mori','Brugia_malayi','Caenorhabditis_brenneri','Caenorhabditis_briggsae','Caenorhabditis_japonica','Caenorhabditis_remanei','Capitella_teleta','Crassostrea_gigas','Culex_quinquefasciatus','Danaus_plexippus','Daphnia_pulex','Dendroctonus_ponderosae','Drosophila_ananassae','Drosophila_erecta','Drosophila_grimshawi','Drosophila_mojavensis','Drosophila_persimilis','Drosophila_pseudoobscura','Drosophila_sechellia','Drosophila_simulans','Drosophila_virilis','Drosophila_willistoni','Drosophila_yakuba','Heliconius_melpomene','Helobdella_robusta','Lepeophtheirus_salmonis','Lucilia_cuprina','Mayetiola_destructor','Melitaea_cinxia','Mnemiopsis_leidyi','Nematostella_vectensis','Octopus_bimaculoides','Pediculus_humanus','Pristionchus_pacificus','Rhodnius_prolixus','Solenopsis_invicta','Stegodyphus_mimosarum','Strigamia_maritima','Tribolium_castaneum','Trichoplax_adhaerens','Zootermopsis_nevadensis','Acropora_digitifera','Aedes_albopictus','Apis_cerana','Apis_dorsata','Apis_florea','Aplysia_californica','Bactrocera_dorsalis','Bactrocera_oleae','Bemisia_tabaci','Bicyclus_anynana','Branchiostoma_belcheri','Branchiostoma_floridae','Camponotus_floridanus','Ceratina_calcarata','Ceratitis_capitata','Ceratosolen_solmsi','Cimex_lectularius','Crassostrea_virginica','Cyphomyrmex_costatus','Drosophila_arizonae','Drosophila_biarmipes','Drosophila_bipectinata','Drosophila_busckii','Drosophila_elegans','Drosophila_eugracilis','Drosophila_ficusphila','Drosophila_kikkawai','Drosophila_miranda','Drosophila_navojoa','Drosophila_obscura','Drosophila_rhopaloa','Drosophila_serrata','Drosophila_suzukii','Drosophila_takahashii','Dufourea_novaeangliae','Echinococcus_granulosus','Eufriesea_mexicana','Eurytemora_affinis','Metaseiulus_occidentalis','Habropoda_laboriosa','Hydra_vulgaris','Limulus_polyphemus','Linepithema_humile','Mizuhopecten_yessoensis','Monomorium_pharaonis','Musca_domestica','Myzus_persicae','Nicrophorus_vespilloides','Nilaparvata_lugens','Onthophagus_taurus','Montastraea_faveolata','Orussus_abietinus','Papilio_machaon','Papilio_polytes','Papilio_xuthus','Parasteatoda_tepidariorum','Pieris_rapae','Pogonomyrmex_barbatus','Priapulus_caudatus','Pseudomyrmex_gracilis','Saccoglossus_kowalevskii','Stomoxys_calcitrans','Trachymyrmex_cornetzi','Trachymyrmex_septentrionalis','Trachymyrmex_zeteki','Varroa_destructor','Vollenhovia_emeryi','Wasmannia_auropunctata','Bactrocera_cucurbitae')
#age cutoffs from James et. al : 2101 for ancient, and 1496 for recent

totalrecent<-total[ which(total$Age_oldest<=1496), ]
totalancient<-total[ which(total$Age_oldest>=2101), ]

totaltransmembrane<-total[which(total$TmhmmTopology >=0.5),]
totalnontransmembrane<- total[which(total$TmhmmTopology <0.5),]



####tranforms

#bc.transform <- function(x,L){(x^L-1)/L}
#asin.transform <- function(x){asin(sqrt(x))}

##test transforms ISD
#ISD.transform <- asin.transform((df$ISD))
#df$ISD.transform <-ISD.transform

#test transforms Clustering
#bc <- boxcox(((df$Clustering)+0.5) ~ 1, lambda = seq(.1,.7,0.01))
#lambda <- bc$x[which.max(bc$y)]

#Clustering.transform <- bc.transform(df$Clustering,lambda)
#df$Clustering.transform<-Clustering.transform

#####subset data

cleanvertebratedata <-total[total$NewickName %in% cleanvertebratelist,]
str(cleanvertebratedata)

cleanvertebratedatarecent <-totalrecent[totalrecent$NewickName %in% cleanvertebratelist,]
str(cleanvertebratedatarecent)

cleanvertebratedataancient <-totalancient[totalancient$NewickName %in% cleanvertebratelist,]
str(cleanvertebratedataancient)

######

cleanvertebratedatatransmembrane <-totaltransmembrane[totaltransmembrane$NewickName %in% cleanvertebratelist,]
str(cleanvertebratedataancient)

cleanvertebratedatanontransmembrane <-totalnontransmembrane[totalnontransmembrane$NewickName %in% cleanvertebratelist,]
str(cleanvertebratedatarecent)


##cleanvertebrate, ISD only

#total pfams
vcISDs1<- lmer(cleanvertebratedata$ISD_WithCys ~ (1|cleanvertebratedata$PfamUID) + cleanvertebratedata$NewickName ,REML=FALSE)
vcsISD1<- summary(vcISDs1)
vcrISD1<- r.squaredGLMM(vcISDs1)
vcaicISD1<- AIC(vcISDs1)

#recent
vcrecISDs1<- lmer(cleanvertebratedatarecent$ISD_WithCys ~ (1|cleanvertebratedatarecent$PfamUID) + cleanvertebratedatarecent$NewickName,REML=FALSE)
vcrecsISD1<- summary(vcrecISDs1)
vcrrecISD1<- r.squaredGLMM(vcrecISDs1)
vcrecaicISD1<- AIC(vcrecISDs1)

#ancient
vcancISDs1<- lmer(cleanvertebratedataancient$ISD_WithCys ~ (1|cleanvertebratedataancient$PfamUID) + cleanvertebratedataancient$NewickName,REML=FALSE)
vcancsISD1<- summary(vcancISDs1)
vcrancISD1<- r.squaredGLMM(vcancISDs1)
vcancaicISD1<- AIC(vcancISDs1)

#transmembrane
vctransISDs1<- lmer(cleanvertebratedatatransmembrane$ISD_WithCys ~ (1|cleanvertebratedatatransmembrane$PfamUID) + cleanvertebratedatatransmembrane$NewickName,REML=FALSE)
vctransISD1<- summary(vctransISDs1)

#nontransmembrane
vcnontransISDs1<- lmer(cleanvertebratedatanontransmembrane$ISD_WithCys ~ (1|cleanvertebratedatanontransmembrane$PfamUID) + cleanvertebratedatanontransmembrane$NewickName,REML=FALSE)
vcnontranssISD1<- summary(vcnontransISDs1)


