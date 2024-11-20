# The purpose of this script is to produce a heatmap illustrating the evolutionary correlations between 
# different components of the avian limb skeleton in altricial and precocial birds. 

library(ggplot2)
 Subplot.title.font.size <- 30/.pt
 Subplot.axis.font.size <- 16/.pt
 Legend.axis.font.size <- 16/.pt
 Legend.title.font.size <- 16/.pt
 Legend.tick.size <- 1
 Subplot.module.linewidth <- 1
 Subplot.tick.size <-1
 Subplot.border.linewidth <- 1
 legend.key.width <- 1
 legend.key.height <- 0.5

setwd('Z:/Andrew_backup/Precocial_paper/all_the_data')
# Set the work directory location to Bjarnason's birds

library( geomorph )
library( ape )
library( nlme )
library(phytools)
# Load required analytical R packages.

load('Csize.22.10.2022.RData')
# Load centroid sizes for study birds. 

load('tree.22.10.2022.RData')
# Load phylogenetic tree of birds
# (Pruned from Prum et al., 2015)
# The object is a phylogeny called 'pruned.tree'.

load('tree.names.22.10.2022.RData')
# Load the tree names required to match birds to the closest genus on the tree.
# The object is a vector called 'tree_names'.

load('masses.22.10.2022.Rdata')
# The object is a named numeric vector called 'masses'.
# The names will not yet match those on the phylogeny.

names(masses) <- tree_names
for( i in 1:length(GPA.Csize) ){
	names(GPA.Csize[[i]]) <- tree_names 
}
# Set names of all properties to be consistent across the sample birds and the phylogeny by matching
# physical specimens to congeners in the phylogeny

get.item <- function( X , item ) { X[[ item ]] }
# This is a simple indexing function

# This is a function that removes allometric scaling/size dependent scaling of limb skeletal lengths from 
# the original data, so that we arrive at a 'deviation from expected' residual.
# How these residuals change with respect to one another tells us about their evolutionary correlation
get.residual.Csize <- function( array, masses, phylogeny, taxa ){ # The function takes a shape array, mass vector, phylogeny and list of taxa
	species<-taxa
	allometry.Csize <- list() # Dumby variable to receive allometric models
	newphy <- drop.tip( phylogeny , phylogeny$tip.label[ !phylogeny$tip.label %in% taxa ] ) # Phylogeny pruned to desired taxa
	for(i in 1:length(array) ){ # For each skeletal element 
		df <- data.frame( mass= log10(masses[taxa]), Csize = log10( array[[ i ]][taxa] ), species = taxa  ) # Define a data frame

			lambda <- phylosig(tree= newphy, x=log10( array[[ i ]][taxa]),  method='lambda')[[1]]
			if(lambda>1){
				lambda<-1
			}
		
		allometry.Csize[[ i ]] <- gls( Csize ~ mass, correlation = corPagel( lambda, phy=newphy, form= ~ species ), data=df ) # Compute an allometric model
	}
	names(allometry.Csize) <- names(array) # Ensure names of allometric models match the skeletal elements they describe 
	residual_Csize <- lapply( lapply( allometry.Csize , get.item , item = "residuals" ) , c ) # Extract model residuals
	return(residual_Csize) # Return residuals
} # Conclude function

matched.birds <- read.csv('precocial_and_altricial_mass_matched_pairs_10_24_2023.csv')
# Load a dataset naming altricial and precocial birds with similar body masses
# It is important we use mass-matched pairs in this study because previous work has shown that variety
# in body mass can influence evolutionary correlations, (https://doi.org/10.1038/s41467-024-48324-y)
# so we do not want body mass variety to be convolved with developmental strategy. 
matched.birds <- matched.birds[-c(28),]
# snip row 28 because it has Spheniscus, which is likely to be a highly influential point.


precocial.tree <- keep.tip(pruned.tree,matched.birds[,1])
altricial.tree <- keep.tip(pruned.tree,matched.birds[,3])
# Phylogenies pruned from Prum et al., 2015 to just the taxa of interest within each.

#names(GPA.Csize) # elements 6,8,9,11,12,13

residual.Csize.precocial <- get.residual.Csize( array = GPA.Csize[c(6,8,9,11,12,13)], masses = masses[matched.birds[,1]], phylogeny = precocial.tree, taxa=matched.birds[,1] )
residual.Csize.altricial <- get.residual.Csize( array = GPA.Csize[c(6,8,9,11,12,13)], masses = masses[matched.birds[,3]], phylogeny = altricial.tree, taxa=matched.birds[,3] )
# Calculate the residual relative sizes of bones once allometric scaling patterns have been removed.

elementsa <- names(GPA.Csize[c(6,8,9,11,12,13)])
pairs <- combn(elementsa,2)
# All possible pairwise combinations of bones, between which evolutionary correlations might be computed

precocial.battenberga<- matrix(NA,length(elementsa),length(elementsa))
rownames(precocial.battenberga)<-colnames(precocial.battenberga)<-elementsa
diff.battenberga <- altricial.battenberga<- precocial.battenberga 
diff.battenberga.p <- altricial.battenberga.p<- precocial.battenberga.p <- precocial.battenberga
# Prepare matrices to recieve results. 

keeps <- list()
ints <-list()

# The following for loop will compute all pairwise combinations of evolutionary correlations in subsets of 
# mass-matched precocial and altricial birds: 

for(i in 1:dim(pairs)[2]){
	bone1 <- pairs[1,i]
	bone2 <- pairs[2,i]
	int <- phylo.integration( residual.Csize.precocial[[bone1]][matched.birds[,1]][precocial.tree$tip] , residual.Csize.precocial[[bone2]][matched.birds[,1]][precocial.tree$tip], phy= precocial.tree)
	
	row <- which(rownames(precocial.battenberga)==bone1)
	column <- which(colnames(precocial.battenberga)==bone2)
	precocial.battenberga[row,column] <- int$Z[[1]]/sqrt(length(matched.birds[,1])) 
	precocial.battenberga.p[row,column] <- int$P.value[[1]]

	row <- which(rownames(precocial.battenberga)==bone2)
	column <- which(colnames(precocial.battenberga)==bone1)
	precocial.battenberga[row,column] <- int$Z[[1]]/sqrt(length(matched.birds[,1]))
	precocial.battenberga.p[row,column] <- int$P.value[[1]]		
	keeps[[i]] <- keep <- int

	ints[[i]] <- int <- phylo.integration( residual.Csize.altricial[[bone1]][matched.birds[,3]][altricial.tree$tip] , residual.Csize.altricial[[bone2]][matched.birds[,3]][altricial.tree$tip], phy= altricial.tree)

	row <- which(rownames(altricial.battenberga)==bone1)
	column <- which(colnames(altricial.battenberga)==bone2)
	altricial.battenberga[row,column] <- int$Z[[1]]/sqrt(length(matched.birds[,3]))
	altricial.battenberga.p[row,column] <- int$P.value[[1]]

	row <- which(rownames(altricial.battenberga)==bone2)
	column <- which(colnames(altricial.battenberga)==bone1)
	altricial.battenberga[row,column] <- int$Z[[1]]/sqrt(length(matched.birds[,3]))
	altricial.battenberga.p[row,column] <- int$P.value[[1]]

	comp <- compare.pls(keep,int)

	row <- which(rownames(diff.battenberga)==bone2)
	column <- which(colnames(diff.battenberga)==bone1)
	diff.battenberga[row,column] <- comp[[3]][1,2]
	diff.battenberga.p[row,column] <- comp[[4]][1,2]

	row <- which(rownames(diff.battenberga)==bone1)
	column <- which(colnames(diff.battenberga)==bone2)
	diff.battenberga[row,column] <- comp[[3]][1,2]
	diff.battenberga.p[row,column] <- comp[[4]][1,2]

	
	print(i/length(pairs))

}


library(reshape2)
# This package contains functions for formatting data 
pdfa <- melt(precocial.battenberga)
pdfpa <- melt(precocial.battenberga.p)
pdfa[which(pdfpa[,3]>0.05),3]<-0
colnames(pdfa)<-c('bone1','bone2','Z')
# Prepare precocial statistical results

adfa <- melt(altricial.battenberga)
adfpa <- melt(altricial.battenberga.p)
adfa[which(adfpa[,3]>0.05),3]<-0
colnames(adfa)<-c('bone1','bone2','Z')
# Prepare altricial statistical results

limits.one <- c(0,1)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# This is a colour blind friendly palette. 

bone_coloursa<-rev(c(rep(cbbPalette[4],3),rep(cbbPalette[8],3)))
# Colours for different modules.



axis<-
ggplot(data=adfa)+
geom_tile(aes(x=NA, y=factor(bone2, levels= elementsa[(c(1:6))])),fill='white')+
scale_y_discrete( limits= rev )+
scale_x_discrete(labels=elementsa[(c(1:6))], guide = guide_axis(angle = 90))+
scale_fill_gradient2(high='black',low='white',midpoint=0.0,na.value='white', limits=c(0,1) )+
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
                     axis.text.y=element_text(size=15,colour=(bone_coloursa),face = "bold"),
				axis.title.x=element_blank(),
				axis.title.y=element_blank(),
				legend.key.width=unit(1/2,'cm'),
				legend.key.height=unit(3,'cm'), 
				legend.title=element_text(size=15),
				legend.text=element_text(size=15),
				plot.margin=unit(c(0.1,0.1,0.1,0.1),'cm'),
                     plot.title=element_text(size=25,hjust=.5), legend.position="bottom",panel.background = element_blank())
# This is a plot element; it is simply an axis with labels. 


precocial.plot.bjarnason<-
ggplot(data=pdfa)+
geom_tile(aes(x=factor(bone1, levels= elementsa), y=factor(bone2, levels= elementsa), fill=Z))+
scale_y_discrete( limits= rev )+
scale_x_discrete(labels=elementsa, guide = guide_axis(angle = 90))+
scale_fill_gradient2(high='black',low='white',midpoint=0.0,na.value='white', limits=limits.one )+
theme(axis.text.x=element_text(size=Subplot.axis.font.size, angle=90, vjust=0.3,colour=rev(bone_coloursa),face = "bold"),
                     #axis.text.y=element_text(size=15,colour=(bone_coloursa),face = "bold"),
				axis.text.y=element_blank(),
				axis.title.x=element_blank(),
 				axis.ticks=element_line(size=2),
				axis.title.y=element_blank(),
				legend.key.width=unit(legend.key.width,'cm'),
				legend.key.height=unit(legend.key.height,'cm'),  
				legend.title=element_text(size=Legend.title.font.size),
				legend.text=element_text(size=Legend.axis.font.size),
				plot.margin=unit(c(0.1,0.1,0.1,0.1),'cm'),
                     plot.title=element_text(size=Subplot.title.font.size,hjust=.5), legend.position="bottom",panel.background = element_blank())+
	labs(fill = expression(paste( italic('Z/'),sqrt(n)  )) )+
geom_segment(data=as.data.frame(cbind(x=c(0.5,0.5,0.5,0.5,3.5,6.5),xend=c(6.5,6.5,6.5,0.5,3.5,6.5),y=c(0.5,3.5,6.5,0.5,0.5,0.5),yend=c(0.5,3.5,6.5,6.5,6.5,6.5) )), 
aes(x=x,xend=xend,y=y,yend=yend), size=1, colour='black')+
annotate('rect',colour='black',fill=NA,xmin=3.5,xmax=6.5,ymin=0.5,ymax=3.5,linewidth=Subplot.module.linewidth+1/4)+
annotate('rect',colour=cbbPalette[8],fill=NA,xmin=3.5,xmax=6.5,ymin=0.5,ymax=3.5,linewidth=Subplot.module.linewidth)+
annotate('rect',colour='black',fill=NA,xmin=0.5,xmax=3.5,ymin=3.5,ymax=6.5,linewidth=Subplot.module.linewidth)+
annotate('rect',colour=cbbPalette[4],fill=NA,xmin=0.5,xmax=3.5,ymin=3.5,ymax=6.5,linewidth=Subplot.module.linewidth+1/4)+
ggtitle('Precocial n=77')+
 guides(fill = guide_colourbar( ticks.linewidth=Legend.tick.size, frame.colour = 'white',
  frame.linewidth = 0/.pt))+ coord_fixed()

battenberg_legend <- cowplot::get_legend(precocial.plot.bjarnason)
# Extract the legend for future plotting. 

elements <- elementsa



precocial.plot.bjarnason<-
ggplot(data=pdfa)+
geom_tile(aes(x=factor(bone1, levels= elements[(c(1:6))]), y=factor(bone2, levels= elements[(c(1:6))]), fill=Z))+
scale_y_discrete( limits= rev,labels=rev(substr(elements[(c(1:6))],1,2)), guide = guide_axis(angle = 0) )+
scale_x_discrete(labels=substr(elements[(c(1:6))],1,2), guide = guide_axis(angle = 90))+
scale_fill_gradient2(high='black',low='white',midpoint=0.0,na.value='white', limits=c(0,1) )+
theme(axis.text.x=element_text(size=Subplot.axis.font.size, angle=90, vjust=0.3,colour=rev(bone_coloursa),face = "bold"),
	axis.text.y=element_text(size=Subplot.axis.font.size, angle=90, vjust=0.3,colour=(bone_coloursa),face = "bold"),
                     #axis.text.y=element_text(size=15,colour=(bone_coloursa),face = "bold"),
				#axis.text.y=element_blank(),
				plot.background = element_rect(fill='transparent', color=NA),
				axis.title.x=element_blank(),
 				axis.ticks=element_line(size=Subplot.tick.size),
				axis.title.y=element_blank(),
				legend.key.width=unit(legend.key.width,'cm'),
				legend.key.height=unit(legend.key.height,'cm'), 
				legend.title=element_text(size=Legend.title.font.size),
				legend.text=element_text(size=Legend.axis.font.size ),
				plot.margin=unit(c(0.1,0.1,0.1,0.1),'cm'),
                     plot.title=element_text(colour= "#56B4E9", size=Subplot.title.font.size,hjust=.5), legend.position="none",panel.background = element_blank())+
	labs(fill = expression(paste( italic('Z/'),sqrt(n)  )) )+
geom_segment(data=as.data.frame(cbind(x=c(0.5,0.5,0.5,0.5,3.5,6.5),xend=c(6.5,6.5,6.5,0.5,3.5,6.5),y=c(0.5,3.5,6.5,0.5,0.5,0.5),yend=c(0.5,3.5,6.5,6.5,6.5,6.5) )), 
aes(x=x,xend=xend,y=y,yend=yend), size=1, colour='black')+
annotate('rect',colour='black',fill=NA,xmin=3.5,xmax=6.5,ymin=0.5,ymax=3.5,linewidth=Subplot.module.linewidth+1/4 )+
annotate('rect',colour=cbbPalette[8],fill=NA,xmin=3.5,xmax=6.5,ymin=0.5,ymax=3.5,linewidth=Subplot.module.linewidth )+
annotate('rect',colour='black',fill=NA,xmin=0.5,xmax=3.5,ymin=3.5,ymax=6.5,linewidth=Subplot.module.linewidth+1/4 )+
annotate('rect',colour=cbbPalette[4],fill=NA,xmin=0.5,xmax=3.5,ymin=3.5,ymax=6.5,linewidth=Subplot.module.linewidth )+
ggtitle('Precocial n=77')+
 guides(fill = guide_colourbar( ticks.linewidth=Legend.tick.size, frame.colour = 'white',
  frame.linewidth = 0/.pt))+ coord_fixed()



altricial.plot.bjarnason<-
ggplot(data=adfa)+
geom_tile(aes(x=factor(bone1, levels= elements[(c(1:6))]), y=factor(bone2, levels= elements[(c(1:6))]), fill=Z))+
scale_y_discrete( limits= rev,labels=rev(substr(elements[(c(1:6))],1,2)), guide = guide_axis(angle = 0) )+
scale_x_discrete(labels=substr(elements[(c(1:6))],1,2), guide = guide_axis(angle = 90))+
scale_fill_gradient2(high='black',low='white',midpoint=0.0,na.value='white', limits=c(0,1) )+
theme(axis.text.x=element_text(size=Subplot.axis.font.size, angle=90, vjust=0.3,colour=rev(bone_coloursa),face = "bold"),
	axis.text.y=element_text(size=Subplot.axis.font.size, angle=90, vjust=0.3,colour=(bone_coloursa),face = "bold"),
                     #axis.text.y=element_text(size=15,colour=(bone_coloursa),face = "bold"),
				#axis.text.y=element_blank(),
				plot.background = element_rect(fill='transparent', color=NA),
				axis.title.x=element_blank(),
 				axis.ticks=element_line(size=Subplot.tick.size),
				axis.title.y=element_blank(),
				legend.key.width=unit(legend.key.width,'cm'),
				legend.key.height=unit(legend.key.height,'cm'), 
				legend.title=element_text(size=Legend.title.font.size),
				legend.text=element_text(size=Legend.axis.font.size ),
				plot.margin=unit(c(0.1,0.1,0.1,0.1),'cm'),
                     plot.title=element_text(colour= '#E69F00', size=Subplot.title.font.size,hjust=.5), legend.position="none",panel.background = element_blank())+
	labs(fill = expression(paste( italic('Z/'),sqrt(n)  )) )+
geom_segment(data=as.data.frame(cbind(x=c(0.5,0.5,0.5,0.5,3.5,6.5),xend=c(6.5,6.5,6.5,0.5,3.5,6.5),y=c(0.5,3.5,6.5,0.5,0.5,0.5),yend=c(0.5,3.5,6.5,6.5,6.5,6.5) )), 
aes(x=x,xend=xend,y=y,yend=yend), size=1, colour='black')+
annotate('rect',colour='black',fill=NA,xmin=3.5,xmax=6.5,ymin=0.5,ymax=3.5,linewidth=Subplot.module.linewidth+1/4 )+
annotate('rect',colour=cbbPalette[8],fill=NA,xmin=3.5,xmax=6.5,ymin=0.5,ymax=3.5,linewidth=Subplot.module.linewidth )+
annotate('rect',colour='black',fill=NA,xmin=0.5,xmax=3.5,ymin=3.5,ymax=6.5,linewidth=Subplot.module.linewidth+1/4 )+
annotate('rect',colour=cbbPalette[4],fill=NA,xmin=0.5,xmax=3.5,ymin=3.5,ymax=6.5,linewidth=Subplot.module.linewidth )+
ggtitle('Altricial n=34')+
 guides(fill = guide_colourbar( ticks.linewidth=Legend.tick.size, frame.colour = 'white',
  frame.linewidth = 0/.pt))+ coord_fixed()



# We are now going to repeat this exercise for species sourced from Brinkworth et al., 2023's dataset of birds

cooney.tree <- read.tree('Prum_merge_hackett_stage2_1K_mcc.tree')
# Load the compiled tree that Brinkworth et al., 2023 use
brinkworth.data <- read.csv('altriciality_scores_brinkworth_11_14_2024.csv')
# Load the Brinkworth dataset, which has been combined with assessments of developmental mode. 

matched.birds <- read.csv('matched_birds_brinkworth_06_03_2024.csv')
# Read the spreadsheet which contains pairs of altricial and precocial bird names from 
# Brinkworth et al., matched by body mass

morphology <- list()
# Define a list to receive morphological variables 

for(i in 1:6){
	morphology[[i]] <- array(data=NA, dim= c(1,1,dim(brinkworth.data)[1]) )
}
# make sure that this list has the correct dimensions to receive all the data

taxa<-intersect(c(matched.birds[,1],matched.birds[,3]), brinkworth.data$JetzTreeNam) 
# Select those taxa which have closely matched masses (altricial and precocial pairs)
# and which also occur in the phylogeny 

element_indices <- c(5,7,8,9,10,11)+1
# These are the columns of the Brinkworth dataset corresponding to the bones of interest 

for(i in 1:6){
	morphology[[i]] <- rep(NA,dim(matched.birds)[1]*2)
	for(j in 1:length(taxa)){
		morphology[[i]][j] <- brinkworth.data[match(taxa[j], brinkworth.data$JetzTreeNam),element_indices[i]]
	}
	names(morphology[[i]]) <- taxa
}
# Populate the morphology object with values

names(morphology) <- c('humerus','radius','carpometacarpus','femur','tibiotarsus','tarsometatarsus')
# Update names

precocial.tree <- keep.tip(cooney.tree,matched.birds[,1])
altricial.tree <- keep.tip(cooney.tree,matched.birds[,3])
# Prune versions of the Jetz tree to only the precocial and altricial mass-matched species from 
# Brinkworth et al., 2023

masses <- brinkworth.data$Body_mass_g[match(taxa, brinkworth.data$JetzTreeNam)]
names(masses) <- taxa
# Extract body masses from Brinkworth dataset

residual.Csize.precocial <- get.residual.Csize( array = morphology, masses = masses[matched.birds[,1]], phylogeny = precocial.tree, taxa=matched.birds[,1] )
residual.Csize.altricial <- get.residual.Csize( array = morphology, masses = masses[matched.birds[,3]], phylogeny = altricial.tree, taxa=matched.birds[,3] )
# Calculate residual bone lengths with allometry removed. 

elements <- names(morphology)
pairs <- combn(elements,2)
# All pairwise combinations defined

precocial.battenberg<- matrix(NA,length(elements),length(elements))
rownames(precocial.battenberg)<-colnames(precocial.battenberg)<-elements
diff.battenberg <- altricial.battenberg<- precocial.battenberg 
diff.battenberg.p <- altricial.battenberg.p<- precocial.battenberg.p <- precocial.battenberg
# Prepare matrices to receive test statistics

keeps <- list()
ints <-list()

# Perform all tests for both precocial and altricial cohorts:

for(i in 1:dim(pairs)[2]){
	bone1 <- pairs[1,i]
	bone2 <- pairs[2,i]
	int <- phylo.integration( residual.Csize.precocial[[bone1]][matched.birds[,1]][precocial.tree$tip] , residual.Csize.precocial[[bone2]][matched.birds[,1]][precocial.tree$tip], phy= precocial.tree)
	
	row <- which(rownames(precocial.battenberg)==bone1)
	column <- which(colnames(precocial.battenberg)==bone2)
	precocial.battenberg[row,column] <- int$Z[[1]]/sqrt(length(matched.birds[,1]))
	precocial.battenberg.p[row,column] <- int$P.value[[1]]

	row <- which(rownames(precocial.battenberg)==bone2)
	column <- which(colnames(precocial.battenberg)==bone1)
	precocial.battenberg[row,column] <- int$Z[[1]]/sqrt(length(matched.birds[,1]))
	precocial.battenberg.p[row,column] <- int$P.value[[1]]
		
	keeps[[i]] <- keep <- int

	ints[[i]] <- int <- phylo.integration( residual.Csize.altricial[[bone1]][matched.birds[,3]][altricial.tree$tip] , residual.Csize.altricial[[bone2]][matched.birds[,3]][altricial.tree$tip], phy= altricial.tree)

	row <- which(rownames(altricial.battenberg)==bone1)
	column <- which(colnames(altricial.battenberg)==bone2)
	altricial.battenberg[row,column] <- int$Z[[1]]/sqrt(length(matched.birds[,3]))
	altricial.battenberg.p[row,column] <- int$P.value[[1]]

	row <- which(rownames(altricial.battenberg)==bone2)
	column <- which(colnames(altricial.battenberg)==bone1)
	altricial.battenberg[row,column] <- int$Z[[1]]/sqrt(length(matched.birds[,3]))
	altricial.battenberg.p[row,column] <- int$P.value[[1]]


	comp <- compare.pls(keep,int)

	row <- which(rownames(diff.battenberg)==bone2)
	column <- which(colnames(diff.battenberg)==bone1)
	diff.battenberg[row,column] <- comp[[3]][1,2]
	diff.battenberg.p[row,column] <- comp[[4]][1,2]

	row <- which(rownames(diff.battenberg)==bone1)
	column <- which(colnames(diff.battenberg)==bone2)
	diff.battenberg[row,column] <- comp[[3]][1,2]
	diff.battenberg.p[row,column] <- comp[[4]][1,2]

	
	print(i/length(pairs))

}

pdf <- melt(precocial.battenberg)
pdfp <- melt(precocial.battenberg.p)
pdf[which(pdfp[,3]>0.05),3]<-0
colnames(pdf)<-c('bone1','bone2','Z')

adf <- melt(altricial.battenberg)
adfp <- melt(altricial.battenberg.p)
adf[which(adfp[,3]>0.05),3]<-0
colnames(adf)<-c('bone1','bone2','Z')

# Prepare precocial and altricial dataframes

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# This is a colour blind friendly palette. 
bone_colours<-rev(c(rep(cbbPalette[4],3),rep(cbbPalette[8],3)))
# Colours for different modules.

precocial.plot.brinkworth<-
ggplot(data=pdf)+
geom_tile(aes(x=factor(bone1, levels= elements[(c(1:6))]), y=factor(bone2, levels= elements[(c(1:6))]), fill=Z))+
scale_y_discrete( limits= rev,labels=rev(substr(elements[(c(1:6))],1,2)), guide = guide_axis(angle = 0) )+
scale_x_discrete(labels=substr(elements[(c(1:6))],1,2), guide = guide_axis(angle = 90))+
scale_fill_gradient2(high='black',low='white',midpoint=0.0,na.value='white', limits=c(0,1) )+
theme(axis.text.x=element_text(size=Subplot.axis.font.size, angle=90, vjust=0.3,colour=rev(bone_coloursa),face = "bold"),
	axis.text.y=element_text(size=Subplot.axis.font.size, angle=90, vjust=0.3,colour=(bone_coloursa),face = "bold"),
                     #axis.text.y=element_text(size=15,colour=(bone_coloursa),face = "bold"),
				#axis.text.y=element_blank(),
				plot.background = element_rect(fill='transparent', color=NA),
				axis.title.x=element_blank(),
 				axis.ticks=element_line(size=Subplot.tick.size),
				axis.title.y=element_blank(),
				legend.key.width=unit(legend.key.width,'cm'),
				legend.key.height=unit(legend.key.height,'cm'), 
				legend.title=element_text(size=Legend.title.font.size),
				legend.text=element_text(size=Legend.axis.font.size ),
				plot.margin=unit(c(0.1,0.1,0.1,0.1),'cm'),
                     plot.title=element_text(colour= "#56B4E9", size=Subplot.title.font.size,hjust=.5), legend.position="none",panel.background = element_blank())+
	labs(fill = expression(paste( italic('Z/'),sqrt(n)  )) )+
geom_segment(data=as.data.frame(cbind(x=c(0.5,0.5,0.5,0.5,3.5,6.5),xend=c(6.5,6.5,6.5,0.5,3.5,6.5),y=c(0.5,3.5,6.5,0.5,0.5,0.5),yend=c(0.5,3.5,6.5,6.5,6.5,6.5) )), 
aes(x=x,xend=xend,y=y,yend=yend), size=1, colour='black')+
annotate('rect',colour='black',fill=NA,xmin=3.5,xmax=6.5,ymin=0.5,ymax=3.5,linewidth=Subplot.module.linewidth+1/4 )+
annotate('rect',colour=cbbPalette[8],fill=NA,xmin=3.5,xmax=6.5,ymin=0.5,ymax=3.5,linewidth=Subplot.module.linewidth )+
annotate('rect',colour='black',fill=NA,xmin=0.5,xmax=3.5,ymin=3.5,ymax=6.5,linewidth=Subplot.module.linewidth+1/4 )+
annotate('rect',colour=cbbPalette[4],fill=NA,xmin=0.5,xmax=3.5,ymin=3.5,ymax=6.5,linewidth=Subplot.module.linewidth )+
ggtitle('Precocial n=77')+
 guides(fill = guide_colourbar( ticks.linewidth=Legend.tick.size, frame.colour = 'white',
  frame.linewidth = 0/.pt))+ coord_fixed()

altricial.plot.brinkworth<-
ggplot(data=adf)+
geom_tile(aes(x=factor(bone1, levels= elements[(c(1:6))]), y=factor(bone2, levels= elements[(c(1:6))]), fill=Z))+
scale_y_discrete( limits= rev,labels=rev(substr(elements[(c(1:6))],1,2)), guide = guide_axis(angle = 0) )+
scale_x_discrete(labels=substr(elements[(c(1:6))],1,2), guide = guide_axis(angle = 90))+
scale_fill_gradient2(high='black',low='white',midpoint=0.0,na.value='white', limits=c(0,1) )+
theme(axis.text.x=element_text(size=Subplot.axis.font.size, angle=90, vjust=0.3,colour=rev(bone_coloursa),face = "bold"),
	axis.text.y=element_text(size=Subplot.axis.font.size, angle=90, vjust=0.3,colour=(bone_coloursa),face = "bold"),
                     #axis.text.y=element_text(size=15,colour=(bone_coloursa),face = "bold"),
				#axis.text.y=element_blank(),
				plot.background = element_rect(fill='transparent', color=NA),
				axis.title.x=element_blank(),
 				axis.ticks=element_line(size=Subplot.tick.size),
				axis.title.y=element_blank(),
				legend.key.width=unit(legend.key.width,'cm'),
				legend.key.height=unit(legend.key.height,'cm'), 
				legend.title=element_text(size=Legend.title.font.size),
				legend.text=element_text(size=Legend.axis.font.size ),
				plot.margin=unit(c(0.1,0.1,0.1,0.1),'cm'),
                     plot.title=element_text(colour= '#E69F00', size=Subplot.title.font.size,hjust=.5), legend.position="none",panel.background = element_blank())+
	labs(fill = expression(paste( italic('Z/'),sqrt(n)  )) )+
geom_segment(data=as.data.frame(cbind(x=c(0.5,0.5,0.5,0.5,3.5,6.5),xend=c(6.5,6.5,6.5,0.5,3.5,6.5),y=c(0.5,3.5,6.5,0.5,0.5,0.5),yend=c(0.5,3.5,6.5,6.5,6.5,6.5) )), 
aes(x=x,xend=xend,y=y,yend=yend), size=1, colour='black')+
annotate('rect',colour='black',fill=NA,xmin=3.5,xmax=6.5,ymin=0.5,ymax=3.5,linewidth=Subplot.module.linewidth+1/4 )+
annotate('rect',colour=cbbPalette[8],fill=NA,xmin=3.5,xmax=6.5,ymin=0.5,ymax=3.5,linewidth=Subplot.module.linewidth )+
annotate('rect',colour='black',fill=NA,xmin=0.5,xmax=3.5,ymin=3.5,ymax=6.5,linewidth=Subplot.module.linewidth+1/4 )+
annotate('rect',colour=cbbPalette[4],fill=NA,xmin=0.5,xmax=3.5,ymin=3.5,ymax=6.5,linewidth=Subplot.module.linewidth )+
ggtitle('Altricial n=77')+
 guides(fill = guide_colourbar( ticks.linewidth=Legend.tick.size, frame.colour = 'white',
  frame.linewidth = 0/.pt))+ coord_fixed()

blank<-ggplot()+ geom_blank()+
  theme(
	axis.text.x=element_text(size=12,colour='black'),
  	axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
# Done


# Now it is time to load illustrated visual elements
library(jpeg)
# Package for processing jpeg image files

precocial.img <- readJPEG("precocial_bird.jpeg")
altricial.img <- readJPEG("altricial_bird.jpeg")
# Illustrations loaded

library(cowplot)
# Useful package for combining images and subplots

  ggdraw() +
  draw_plot(blank) +
  draw_image(altricial.img,x=0.37,y=0.10,width=0.8,height=.85)+
  draw_image(precocial.img,x=-0.245,y=0.055,width=0.8,height=.85)+
  draw_plot(precocial.plot.brinkworth, x = 0.265, y = -0.28, width = .2, height = 1)+
  draw_plot(altricial.plot.brinkworth, x = 0.46, y = -0.28, width = .2, height = 1)+ #-0.175
  draw_plot(precocial.plot.bjarnason, x = 0.265, y = 0.08, width = .2, height = 1)+
  draw_plot(altricial.plot.bjarnason, x = 0.46, y = 0.08, width = .2, height = 1)+ # 0.235
	draw_plot(battenberg_legend,x=0.375,y=0.7,width=0.2,height=0.35)
 # Total compiled plot

#setwd('C:/Users/Lab/Documents/AOrkney/Precocial_study')
setwd('Z:/Andrew_backup/Precocial_paper') 
ggsave(filename='Battenbergs_10_28_2024.pdf',width=18,height=11,unit='cm') 
# save



######
# SCRAPS: 

library(ggpubr)

ggarrange(precocial.plot.bjarnason,altricial.plot.bjarnason,precocial.plot.brinkworth,altricial.plot.brinkworth,align='h',ncol=2,nrow=2, 
labels=c('a','b','c','d'), font.label = list(size = 30, color = "black", face = "bold", family = NULL), hjust=-0.4,common.legend=T )

setwd('C:/Users/Lab/Documents/AOrkney/Precocial_study')
ggsave( height=30,width=45,units='cm',dpi=300,filename='Figure_I_06_04_2024.pdf')


"Fregata_magnificens" "Puffinus_griseus"    "Phaethon_rubricauda"
[4] "Arenaria_interpres"  "Rollandia_rolland"   "Aythya_valisineria" 

# There are only 6 animals that occur in both datasets 

# you got as far as.. Bambusicola_fytchii
