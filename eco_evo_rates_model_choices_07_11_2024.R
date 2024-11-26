# The purpose of this script is to determine whether altricial developmental strategy within the Bjarnason dataset of birds
# facilitates more rapid evolutionary divergence of wing and leg sizes and hence access to novel flight-style and foot-use ecologies
# that were not recorded in precocial birds within the Bjarnason dataset. 
# 'Novelty' of flight-style and foot-use ecology will be mapped as a stochastic character over the phylogeny of Prum et al., 2015
# for species that possess an altricial developmental strategy.
# Thereafter, we will seek to determine whether single-rate Brownian or multiple-rate Brownian motion models depending on evolutionary novelty
# are preferred explanations for the diversity of wing to leg size ratios observed across altricial birds. 

library(geiger) # 2.0.10 Fit 
library( ape ) # 5.7-1
library(phytools) # 1.5-1
library( geomorph ) # 4.0.5
library(mvMORPH) # 1.1.8

setwd('C:/Users/Lab/Documents/AOrkney/bird_scripts_and_data_15_02_2023')
metadata <- read.csv('altriciality_scores_10_15_2023.2.csv')
# Load bird metadata containing developmental strategy classifications and congener matches to Prum phylogeny

# Load some requisite datasets
load('tree.22.10.2022.RData')
load('tree.names.22.10.2022.RData')
load('Csize.22.10.2022.RData')
load('masses.22.10.2022.RData')
# Done

names(masses)<-tree_names
# Ensure names match the phylogeny. 
for( i in 1:length(GPA.Csize) ){
	names(GPA.Csize[[i]]) <- tree_names 
}
# Done
Csize.birds<-GPA.Csize
# Store data in a new object 

developmental_mode <- metadata$my[match(tree_names,metadata$tree)]
developmental_mode[grep('emi', developmental_mode)]<-'Semi'
developmental_mode[grep('tricial', developmental_mode)]<-'Altricial'
developmental_mode[grep('cocial', developmental_mode)]<-'Precocial'
names(developmental_mode)<-tree_names
# Make developmental categorisations consistent and match their names to the Prum phylogeny

tips <- names(masses)[which(developmental_mode == 'Precocial' | developmental_mode == 'Altricial')]
short.tree <- keep.tip(pruned.tree,tips)
short.mode <- developmental_mode[short.tree$tip]
# Restrict developmental modes classifications to only hose taxa we can confidently identify as either precocial or altricial 

leg <- (GPA.Csize$tibiotarsus+GPA.Csize$femur+GPA.Csize$tarsometatarsus)
wing <- (GPA.Csize$humerus+GPA.Csize$radius+GPA.Csize$carpometacarpus)
total <- wing+leg
trait <- ((wing-leg)/total)
# Define the trait as the difference between the wing and leg, as compared to their total size
# This is a measure of limb skeleton 'evenness'




setwd('C:/Users/Lab/Documents/AOrkney/bird_scripts_and_data_15_02_2023')
# Set work directory to location wherein ecological metadata is stored 

Foot.use <- read.csv('standardised_foot_scores_11_07_2023.csv')
# Load Brigit Tronrud's foot-use dataset 

Flight.style <- read.csv('flight_masses_22_10_2022_plus_A.csv')
# Load the flight-style categorisations 


# This loop will identify which foot-uses are not present in any precocial species in the dataset 
for(i in 5:17){
	if(
	length(grep('Precocial',developmental_mode[Foot.use$X[which(Foot.use[i]==1)]]))==0
	){
		print(i)
	}
}
colnames(Foot.use[9:13])
# Hannging_Clinging, Scansorial, Raptorial, Digit_manip and Force_grip occur only in 
# Altricial birds

Novel_feet <- Foot.use$X[which(rowSums(Foot.use[,9:13])>0)]
# These are the birds that practice novel foot ecologies that are not available to precocial birds

# Undertake the same task for flight-styles 
for(i in 3:13){
	if(
	length(grep('Precocial',developmental_mode[Flight.style$X[which(Flight.style[i]==1)]]))==0
	){
		print(i)
	}
}
colnames(Flight.style[c(3,6,11,13)])
# Exclusive_pelagic_soaring, Sallying_flight, Flap_bounding and Hovering are novel ecologies
# that are not available to precocial birds

Novel_flight <- Flight.style$X[which(rowSums(Flight.style[,c(3,6,11,13)])>0)]

altricial.tree <- keep.tip(short.tree,names(developmental_mode)[which(developmental_mode=='Altricial')] )
# Define a phylogeny for only species practicing altricial developmental strategies 

Novel_flight <- Novel_flight[Novel_flight%in%altricial.tree$tip]
Novel_feet <- Novel_feet[Novel_feet%in%altricial.tree$tip]
# Ensure that taxa practicing novel foot-uses and flight-styles occur in the altricial tree 

altricial.mode <- rep('conservative',length(altricial.tree$tip))
names(altricial.mode)<-altricial.tree$tip
altricial.mode[Novel_flight]<-'novel'
# Codify ecological 'novelty' of flight-style across the tips of the phylogeny relating altricial birds. 

dat <- rep(1,length(altricial.tree$tip))
names(dat) <- altricial.tree$tip
dat[Novel_flight]<-2
# Render the states as numbers

#model = c("ER","SYM","ARD","meristic") # I do not believe we need this line 

novel.flight.ER <- fitDiscrete(phy=altricial.tree, dat=dat, model='ER')
novel.flight.ER$opt$aic
novel.flight.SYM <- fitDiscrete(phy=altricial.tree, dat=dat, model='SYM')
novel.flight.SYM$opt$aic
novel.flight.ARD <- fitDiscrete(phy=altricial.tree, dat=dat, model='ARD')
novel.flight.ARD$opt$aic
novel.flight.meristic <- fitDiscrete(phy=altricial.tree, dat=dat, model='meristic')
# Problem with meristic model fit?
novel.flight.meristic$opt$aic
# The ARD model has the preferred AIC, but none of the AIC diffs are greater than 2 anyway

alt.tree<-make.simmap(altricial.tree,altricial.mode[which(names(altricial.mode) %in% altricial.tree$tip)] , model="ARD", nsim=1, pi= c(0,1) )
alt.tree$log
alt.tree<-make.simmap(altricial.tree,altricial.mode[which(names(altricial.mode) %in% altricial.tree$tip)] , model="ARD", nsim=1, pi= c(1,0) )
alt.tree$log
# The log likelihood is higher when novel flight styles are assumed to have originated ancestrally in altricial birds. 

alt.tree<-make.simmap(altricial.tree,altricial.mode[which(names(altricial.mode) %in% altricial.tree$tip)] , model="ARD", nsim=1000,  pi= c(0,1) )
col<-c("orange","black"); names(col)<-c("conservative","novel")
plotSimmap(alt.tree[[1]],col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)
# An All-Rates-Different stochastic character map of 'novelty'
# with novel flight-styles assumed ancestral to altricials, produces a stable reconstruction

# We will now determine whether the inferred rate of divergent wing and leg evolution across the altricial phylogeny is 
# explained by flight-style novelty, by comparing single and multiple-rate Brownian motion models. 

preferred.model<-list()
conservative.rate<-list()
novel.rate<-list()
AIC.dif <- list()
for(i in 1:length(alt.tree)){
	model1.alt <- mvBM(tree=alt.tree[[i]], data=trait[alt.tree[[i]]$tip], model='BM1', optimization ="subplex")
	model2.alt <- mvBM(tree=alt.tree[[i]], data= trait[alt.tree[[i]]$tip] , model='BMM', optimization ="subplex")
	if(model1.alt$conv==0 & model1.alt$hess.value==0 & model2.alt$conv==0 & model2.alt$hess==0){
		conservative.rate[[i]]<-model2.alt$sigma[[1]]
		novel.rate[[i]]<-model2.alt$sigma[[2]]
		results <- list(model1.alt,model2.alt)
		weights <- mvMORPH::aicw(results)
		preferred.model[[i]] <- which(weights$diff==0)
		AIC.dif[[i]] <- max(weights$diff)
		print(paste((i*100)/length(alt.tree),'% done'))
	}
}

table(unlist(preferred.model))
# There is a universal preference for a two-rate model 
median(unlist(conservative.rate))/median(unlist(novel.rate))
# The altricials exhibitting novel ecologies evolve much more rapidly. 
plot(unlist(conservative.rate),unlist(novel.rate))
# There doesn't appear to be any relationship between the inferred rate magnitudes. 
range(unlist(AIC.dif))
# delta AIC is always high 

# Summary statistics: 
median(unlist(novel.rate)[which(unlist(AIC.dif)>=2)])
mean(unlist(novel.rate)[which(unlist(AIC.dif)>=2)])
sd(unlist(novel.rate)[which(unlist(AIC.dif)>=2)])

median(unlist(conservative.rate)[which(unlist(AIC.dif)>=2)])
mean(unlist(conservative.rate)[which(unlist(AIC.dif)>=2)])
sd(unlist(conservative.rate)[which(unlist(AIC.dif)>=2)])

# The same approach will now be repeated for novel foot-use ecologies: 

dat <- rep(1,length(altricial.tree$tip))
names(dat) <- altricial.tree$tip
dat[Novel_feet]<-2

model = c("ER","SYM","ARD","meristic")

novel.feet.ER <- fitDiscrete(phy=altricial.tree, dat=dat, model='ER')
novel.feet.ER$opt$aic
novel.feet.SYM <- fitDiscrete(phy=altricial.tree, dat=dat, model='SYM')
novel.feet.SYM$opt$aic
novel.feet.ARD <- fitDiscrete(phy=altricial.tree, dat=dat, model='ARD')
novel.feet.ARD$opt$aic
novel.feet.meristic <- fitDiscrete(phy=altricial.tree, dat=dat, model='meristic')
novel.feet.meristic$opt$aic
# The ER model has the preferred AIC, but none of the AIC diffs are greater than 2 anyway

altricial.mode <- rep('conservative',length(altricial.tree$tip))
names(altricial.mode)<-altricial.tree$tip
altricial.mode[Novel_feet]<-'novel'

alt.tree<-make.simmap(altricial.tree,altricial.mode[which(names(altricial.mode) %in% altricial.tree$tip)] , model="ER", nsim=1, pi= c(0,1) )
alt.tree$log
alt.tree<-make.simmap(altricial.tree,altricial.mode[which(names(altricial.mode) %in% altricial.tree$tip)] , model="ER", nsim=1, pi= c(1,0) )
alt.tree$log
# The log likelihood is higher when novel foot uses are assumed to be derived within altricial birds. 

alt.tree<-make.simmap(altricial.tree,altricial.mode[which(names(altricial.mode) %in% altricial.tree$tip)] , model="ER", nsim=1000, pi= c(1,0) )
col<-c("orange","black"); names(col)<-c("conservative","novel")
plotSimmap(alt.tree[[1]],col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)
# An Equal Rates model with novel feet being derived within altricials, produces a stable reconstruction

preferred.model<-list()
conservative.rate<-list()
novel.rate<-list()
AIC.dif <- list()
for(i in 1:length(alt.tree)){
	model1.alt <- mvBM(tree=alt.tree[[i]], data=trait[alt.tree[[i]]$tip], model='BM1', optimization ="subplex")
	model2.alt <- mvBM(tree=alt.tree[[i]], data= trait[alt.tree[[i]]$tip] , model='BMM', optimization ="subplex")
	if(model1.alt$conv==0 & model1.alt$hess.value==0 & model2.alt$conv==0 & model2.alt$hess==0){
		conservative.rate[[i]]<-model2.alt$sigma[[1]]
		novel.rate[[i]]<-model2.alt$sigma[[2]]
		results <- list(model1.alt,model2.alt)
		weights <- mvMORPH::aicw(results)
		preferred.model[[i]] <- which(weights$diff==0)
		AIC.dif[[i]] <- max(weights$diff)
		print(paste((i*100)/length(alt.tree),'% done'))
	}
}
table(unlist(preferred.model))
# There is a substantial preference for a two-rate model 
median(unlist(conservative.rate))/median(unlist(novel.rate))
# The altricials exhibitting novel ecologies evolve much more slowly. 
plot(unlist(conservative.rate),unlist(novel.rate))
# There doesn't appear to be any conclusive relationship between the inferred rate magnitudes. 
range(unlist(AIC.dif))
# delta AIC is usually high 

# Summary statistics: 
median(unlist(novel.rate)[which(unlist(AIC.dif)>=2)])
mean(unlist(novel.rate)[which(unlist(AIC.dif)>=2)])
sd(unlist(novel.rate)[which(unlist(AIC.dif)>=2)])

median(unlist(conservative.rate)[which(unlist(AIC.dif)>=2)])
mean(unlist(conservative.rate)[which(unlist(AIC.dif)>=2)])
sd(unlist(conservative.rate)[which(unlist(AIC.dif)>=2)])



