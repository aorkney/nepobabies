library( ape ) # 5.7-1
library(phytools) # 1.5-1
library( geomorph ) # 4.0.5
library(mvMORPH) # 1.1.8

# We will first explore whether it should be possible to identify 
# differences in evolutionary rate as a function of developmental mode in Bjarnason's bird dataset. 
# We already possess evidence which indicates that the wing and leg proportions evolve more independently
# of one another in altricial rather than precocial birds. 
# We might therefore expect that the difference between wing and leg lengths (normalised to total limb lengths)
# should evolve more rapidly in altricial birds and more slowly in precocial birds. 

# We will use a stochastic character mapping approach to map the possible distribution of developmental
# modes across avian phylogeny- generating 1000 possible distributions. 
# We will employ the constraint that the ancestor of the radiation of crown birds is probably precocial.

# We will compare single Brownian and multiple Brownian rate models.
# We will use the delta AIC statistic as the comparison statistic.

# We will also investigate Brinkworth's much larger dataset of birds. 

setwd('Z:/Andrew_backup/Precocial_paper/all_the_data')
metadata <- read.csv('altriciality_scores_10_15_2023.2.csv')
# Load metadata for Bjarnason's birds

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
# Store data in object


developmental_mode <- metadata$my[match(tree_names,metadata$tree)]
developmental_mode[grep('emi', developmental_mode)]<-'Semi'
developmental_mode[grep('tricial', developmental_mode)]<-'Altricial'
developmental_mode[grep('cocial', developmental_mode)]<-'Precocial'
names(developmental_mode)<-tree_names
# Organize developmental strategy classifications

tips <- names(masses)[which(developmental_mode == 'Precocial' | developmental_mode == 'Altricial')]
short.tree <- keep.tip(pruned.tree,tips)
short.mode <- developmental_mode[short.tree$tip]
# Prune the data to the taxa of interest
tree<-make.simmap(short.tree,short.mode , model="ER", nsim=1000, pi= c(0,1) )
# We shall assume an equal rat es model of developmental strategy evolution across the tree, conditioned with a 
# precocial root state.
# The user may change order of pi to assume a different root state
col<-c("blue","orange"); names(col)<-c("Precocial","Altricial")
plotSimmap(tree[[1]],col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)
# The distribution of precocial birds across the tree is compatible with multiple independent origins of altriciality
# different runs of the 'make.simmap' function yield different estimates of the the distribution of developmental strategies through
# evolutionary history. 
# I have set pi as a vector which specifies a 100% certainty for the root state being precocial. 
# I have generated 1000 sample possible ancestral state reconstructions


leg <- (GPA.Csize$tibiotarsus+GPA.Csize$femur+GPA.Csize$tarsometatarsus)
wing <- (GPA.Csize$humerus+GPA.Csize$radius+GPA.Csize$carpometacarpus)
total <- wing+leg
trait <- ((wing-leg)/total)
# Define the trait as the difference between the wing and leg, as compared to their total size
# This is a measure of limb skeleton 'evenness'


preferred.model<-list()
altricial.rate<-list()
precocial.rate<-list()
weights <- list()
for(i in 1:length(tree)){ # For all 1000 possible histories of developmental strategy evolution
	model1 <- mvBM(tree=tree[[i]], data=trait[tree[[i]]$tip], model='BM1', optimization ="subplex") # Compute a single rate brownian model
	model2 <- mvBM(tree=tree[[i]], data= trait[tree[[i]]$tip] , model='BMM', optimization ="subplex") # Compute a multiple rate model depending on development
	if(model1$conv==0 & model1$hess.value==0 & model2$conv==0 & model2$hess==0){
		altricial.rate[[i]]<-model2$sigma[[1]]
		precocial.rate[[i]]<-model2$sigma[[2]]
		results <- list(model1,model2)
		weights[[i]] <- aicw(results) # What are the relative AIC values for the models?
		preferred.model[[i]] <- which(weights[[i]]$diff==0) # Which model has the greatest explanatory power relative to cost?
	}
	print(i)
}
table(unlist(preferred.model))
# A 2-rate model is preferred in Bjarnason's birds, depending on developmental mode 
# across all the sample trees
# It doesn't actually matter if an altricial or precocial root state is assumed 
median(unlist(precocial.rate))/median(unlist(altricial.rate))
# This provides an impression of the relative difference in rates of divergent wing and leg size evolution between precocial and altricial birds

check.diff <- function(input){
	return(abs(input$diff[1]-input$diff[2]))
}
length(which(unlist(lapply(weights, check.diff))>2))/1000
# What proportion of models are preferred with an AIC difference in excess of 2? 

#Summary statistics
mean((unlist(altricial.rate)/unlist(precocial.rate))[which(unlist(lapply(weights, check.diff))>2)])
sd((unlist(altricial.rate)/unlist(precocial.rate))[which(unlist(lapply(weights, check.diff))>2)])
hist((unlist(altricial.rate)/unlist(precocial.rate))[which(unlist(lapply(weights, check.diff))>2)])

mean((unlist(precocial.rate))[which(unlist(lapply(weights, check.diff))>2)])
mean((unlist(altricial.rate))[which(unlist(lapply(weights, check.diff))>2)])

median((unlist(precocial.rate))[which(unlist(lapply(weights, check.diff))>2)])
median((unlist(altricial.rate))[which(unlist(lapply(weights, check.diff))>2)])

sd((unlist(altricial.rate))[which(unlist(lapply(weights, check.diff))>2)])
sd((unlist(precocial.rate))[which(unlist(lapply(weights, check.diff))>2)])

mean((unlist(precocial.rate)))
mean((unlist(altricial.rate)))

sd((unlist(altricial.rate)))
sd((unlist(precocial.rate)))

# 5e-4 +/- 6e-6 rate for altricials
# 2.2e-4 +/- 6e-6 rate for precocials


# Load Brinkworth's aggregated dataset. 
brinkworth.data <- read.csv('altriciality_scores_brinkworth_11_14_2024.csv')

cooney.tree <- read.tree('Prum_merge_hackett_stage2_1K_mcc.tree')
# Load a tree with enough taxa to accommodate Brinkworth's aggregated suite of data

pruned.tree <- keep.tip(cooney.tree, brinkworth.data$JetzTreeName)
# Exclude taxa that are not present in Brinkworth

GPA.Csize <- list()
GPA.Csize[[1]] <- brinkworth.data[,5+1]
GPA.Csize[[2]] <- brinkworth.data[,7+1]
GPA.Csize[[3]] <- brinkworth.data[,8+1]
GPA.Csize[[4]] <- brinkworth.data[,9+1]
GPA.Csize[[5]] <- brinkworth.data[,10+1]
GPA.Csize[[6]] <- brinkworth.data[,11+1]

for(i in 1:length(GPA.Csize)){
	names(GPA.Csize[[i]]) <- brinkworth.data$Jetz
}
# Compile an object of skeletal length measurements from Brinkworth

names(GPA.Csize) <- c('humerus','radius','carpometacarpus','femur','tibiotarsus','tarsometatarsus')
masses <- brinkworth.data$Body
names(masses)<-brinkworth.data$Jetz
# Render names consistent 

metadata <-brinkworth.data
developmental_mode <- metadata$altr[match(names(masses),metadata$Jetz)]
developmental_mode[grep('emi', developmental_mode)]<-'Semi'
developmental_mode[grep('tricial', developmental_mode)]<-'Altricial'
developmental_mode[grep('cocial', developmental_mode)]<-'Precocial'
names(developmental_mode)<-names(masses)
# Compile and make consistent the developmental strategy data 

tips <- names(masses)[which(developmental_mode == 'Precocial' | developmental_mode == 'Altricial')]
short.tree <- keep.tip(pruned.tree,tips)
short.mode <- developmental_mode[short.tree$tip]
# Prune the Brinkworth dataset to only those taxa with clear developmental strategies

tree<-make.simmap(short.tree,short.mode , model="ER", nsim=1000, pi= c(0,1) )
# The user may change the order of pi to change the root state
col<-c("blue","orange"); names(col)<-c("Precocial","Altricial")
plotSimmap(tree[[1]],col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)
# The distribution of precocial birds across the tree is compatible with multiple independent origins of altriciality
# different runs of the 'make.simmap' function yield different estimates of the the distribution of precociality through
# evolutionary history. 
# I have set pi as a vector which specifies a 100% certainty for the root state being precocial. 
# I have generated 1000 possible ancestral state reconstructions

leg <- (GPA.Csize$tibiotarsus+GPA.Csize$femur+GPA.Csize$tarsometatarsus)
wing <- (GPA.Csize$humerus+GPA.Csize$radius+GPA.Csize$carpometacarpus)
total <- wing+leg
trait <- ((wing-leg)/total)
# Define the trait as the difference between the wing and leg, as compared to their total length
# This is a measure of limb skeleton 'evenness'


# Compute evolutionary models of trait over all 1000 possible histories of developmental strategy evolution 
preferred.model<-list()
altricial.rate<-list()
precocial.rate<-list()
weights <- list()
for(i in 1:length(tree)){
	model1 <- mvBM(tree=tree[[i]], data=trait[tree[[i]]$tip], model='BM1', optimization ="subplex") # Compute a single rate brownian model
	model2 <- mvBM(tree=tree[[i]], data= trait[tree[[i]]$tip] , model='BMM', optimization ="subplex") # Compute a multiple rate model depending on development
	if(model1$conv==0 & model1$hess.value==0 & model2$conv==0 & model2$hess==0){
		altricial.rate[[i]]<-model2$sigma[[1]]
		precocial.rate[[i]]<-model2$sigma[[2]]
		results <- list(model1,model2)
		weights[[i]] <- aicw(results) # What are the relative AIC values for the models?
		preferred.model[[i]] <- which(weights[[i]]$diff==0) # Which model has the greatest explanatory power relative to cost?
	}
	print(i)
}
# This loop may take considerable time to run.

length(which(unlist(lapply(weights, check.diff))>2))/1000
# What proportion of the 1000 trees are compatible with a 2-rate model? 


#Summary statistics
mean((unlist(altricial.rate)/unlist(precocial.rate))[which(unlist(lapply(weights, check.diff))>2)])
sd((unlist(altricial.rate)/unlist(precocial.rate))[which(unlist(lapply(weights, check.diff))>2)])
hist((unlist(altricial.rate)/unlist(precocial.rate))[which(unlist(lapply(weights, check.diff))>2)])

mean((unlist(precocial.rate))[which(unlist(lapply(weights, check.diff))>2)])
mean((unlist(altricial.rate))[which(unlist(lapply(weights, check.diff))>2)])

median((unlist(precocial.rate))[which(unlist(lapply(weights, check.diff))>2)])
median((unlist(altricial.rate))[which(unlist(lapply(weights, check.diff))>2)])

sd((unlist(altricial.rate))[which(unlist(lapply(weights, check.diff))>2)])
sd((unlist(precocial.rate))[which(unlist(lapply(weights, check.diff))>2)])

mean((unlist(precocial.rate)))
mean((unlist(altricial.rate)))

sd((unlist(altricial.rate)))
sd((unlist(precocial.rate)))

# ~6.1e-4 +/- 1e-6 rate for altricials
# 3.4e-4 +/- 1e-6 rate for precocials
