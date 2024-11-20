# This script interacts with 
# The supplementary dataset of Brinkworth et al., 2023 (URL)
# It prepares and curates the data for analysis 
# in the study of Orkney et al., 2025

setwd('Z:/Andrew_backup/Precocial_paper/all_the_data')
# Set the work directory to the location of the Brinkworth supplementary download

library(ape) # 5.7-1
# This package contains functions that are used to process phylogenetic trees

cooney.tree <- read.tree('Prum_merge_hackett_stage2_1K_mcc.tree')
# This is a merged tree that was first produced by Cooney, binding dense phylogenies onto the robust
# backbone of the Prum et al., 2015 phylogeny 

brinkworth.data <- read.csv('Supplementary_Data_1.csv')
# Read the limb length measurement data 

dense <- brinkworth.data[complete.cases(brinkworth.data[5:10]),]
# Eliminate birds with incomplete skeletal proportion measurements


# I have noticed that there are occasionally multiple assignments of a species to a single
# JetzTreeName
# What are we to do about these?
# I think we should eliminate the second bird which matches them

for(i in 1:length(unique(dense$JetzTreeName)) ){
	matching.species <- which(dense$JetzTreeName == unique(dense$JetzTreeName)[i]) 
	if(length(matching.species)>1){
		dense <- dense[-matching.species[-1],]
	}
}
# We have eliminated double-matches

families <- dense$Family
unique(families)[order(unique(families))]
# Review the families that are present in Brinkworth's aggregated dataset
# Developmental strategy is conserved at the family level, 
# so we can use familial identity to classify Brinkworth's birds into distinct developmental strategies
# in accordance with Starck et al., URL

altriciality <- rep('NA',length(families))
#Precocial 1:
to_match <- c('Apterygidae', 'Anatidae', 'Anhimidae', 'Jacanidae',
 'Phalaropodidae', 'Pterocledidae',
'Recurvirostridae', 'Charadriidae', 'Scolopacidae',
 'Tinamidae', 'Casuariidae')
# These families are very highly precocial 

altriciality[grep(paste(to_match,collapse='|'),families)] <- 'Precocial_1'

#Precocial 2:
to_match <-c('Rheidae', 'Numididae', 'Phasianidae', 'Tetraonidae',
 'Meleagrididae', 'Struthionidae',
'Dromadidae', 'Rostratulidae', 'Mesoenatidae',
'Aramidae', 'Thinocoridae', 'Heliornithidae',
'Chioniidae')
# These families are highly precocial 

altriciality[grep(paste(to_match,collapse='|'),families)] <- 'Precocial_2'

#Precocial 3:
to_match <-c('Cracidae', 'Gaviidae', 'Glareolidae', 'Gruidae', 'Otidae', 'Psophiidae', 'Rallidae',
'Thrnicidae', 'Pedionomidae', 'Podicipedidae', 'Haematopodidae')
# These families are precocial 

altriciality[grep(paste(to_match,collapse='|'),families)] <- 'Precocial_3'

#Semiprecocial:

to_match <-c('Alcidae', 'Burhinidae', 'Eurypigidae', 'Hydrobatidae',
 'Laridae', 'Opisthocomidae', 'Phoenicopteridae',
 'Rhynchopidae', 'Rhynochetidae', 'Stercorariidae')
# These families are more precocial than not, but have elements of altricial development 

altriciality[grep(paste(to_match,collapse='|'),families)] <- 'Semiprecocial'


to_match <-c('Accipitridae', 'Aegothelidae', 'Ardeidae',
'Balaenicipitidae', 'Caprimulgidae', 'Cariamidae', 
'Cathartidae', 'Falconidae', 'Ciconiidae',
'Dromadidae', 'Hemiprocnidae', 'Musophagidae',
 'Nydibiidae', 'Pandionidae', 'Pelecanoididae',
 'Podargidae', 'Sagittariidae', 'Scopidae',
'Threskiornithidae')
# These families are more altricial than not, but have elements of precocial development

altriciality[grep(paste(to_match,collapse='|'),families)] <- 'Semialtricial'


#Altricial 1:

to_match <-c('Aramidae', 'Columbidae', 'Diomedeidae',
'Phaethontidae', 'Procellariidae', 'Spheniscidae', 'Strigidae', 'Tytonidae')
# These families are highly altricial 

altriciality[grep(paste(to_match,collapse='|'),families)] <- 'Altricial_1'

#Altricial 2:

to_match <-c('Alcedinidae', 'Anhingidae', 'Apodidae', 'Bucconidae',
 'Bucerotidae', 'Capitonidae', 'Coliidae', 'Coraciidae', 'Cuculidae', 'Fregatidae', 
'Galbulidae', 'Indicatoridae', 'Jyngidae',
'Leptosomatidae', 'Meropidae', 'Momotidae',
'Phalacrocoracidae', 'Picidae', 'Psittacidae', 'Rhamphastidae',
 'Sulidae', 'Steatornithidae', 'Todidae', 'Trochilidae', 'Trogonidae', 'Upupidae')
# These families are also unambiguously altricial 

altriciality[grep(paste(to_match,collapse='|'),families)] <- 'Altricial_2'

# Finally, all remaining passeriformes are Altricial 2:

altriciality[grep('Passeriformes',dense$Order)]<-'Altricial_2'

newdata<- cbind(dense,altriciality)
write.csv(x=newdata,file='altriciality_scores_brinkworth_11_14_2024.csv')
# Write the data to a new csv file 