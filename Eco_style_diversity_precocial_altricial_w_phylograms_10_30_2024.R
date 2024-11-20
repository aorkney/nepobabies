# The goal of this script is to:
# Produce phylograms of flight style and foot use variety across the sampled birds
# in Bjarnason's mass-matched dataset...
# Produce ribbon plots for richness accumulation curves showing how the breadth of sampled
# flight-style and foot-use variety increases as a function of sample-size, subset by developmental mode
# Given that both altricials and precocials are widely disseminated through avian phylogeny 
# we can relax concerns that phylogenetic structure drives results. 

library( geomorph ) # 4.0.5
library( ape ) # 5.7-1
library( nlme ) # 3.1-162
library( ggplot2 ) # 3.4.1
library( ggdendro ) # 0.1.23
library( dendextend ) # 1.17.1
library( phytools ) # 1.5-1
library(ggpubr)
library(grid)
# Load required packages

setwd('Z:/Andrew_backup/Precocial_paper/all_the_data')
Bjarnason <- read.csv('precocial_and_altricial_mass_matched_pairs_10_24_2023.csv')
Bjarnason <- Bjarnason[-28,]
# Load the Bjarnason paired birds and remove the pair that includes the penguin. 

load('tree.22.10.2022.RData')
# Load the pruned Bjarnason tree

pruned.tree <- keep.tip(pruned.tree, c(Bjarnason$precocials,Bjarnason$altricials)) 
# Prune the tree further to just the taxa of interest

flight.styles <- read.csv('flight_masses_22_10_2022_plus_A.csv')
footuse.styles <- read.csv('standardised_foot_scores_11_07_2023.csv')
# Great thanks are extended to Brigit Tronrud for compiling this dataset in 2020.
# It is first published in Orkney et al., 2021: URL 
footuse.styles<-cbind(footuse.styles,rep(NA,dim(footuse.styles)[1])) 
# Read the flight and foot-use styles of Bjarnason's birds

# Define a custom function that samples the precocial and altricial populations at
# different sample sizes and totals perceived richness 
sample.rarify<-function(taxa,ecologies,resample.size){

		precocial.variety <- matrix(NA, resample.size*dim(taxa)[1],2)
		precocial.variety[,1] <- rep(1:dim(taxa)[1],each=resample.size)

		altricial.variety <- matrix(NA, resample.size*dim(taxa)[1],2)
		altricial.variety[,1] <- rep(1:dim(taxa)[1],each=resample.size)

	for(i in 1:dim(taxa)[1]){
		for(j in 1:resample.size){
			precocial.subsample<-sample(taxa$precocials,i, replace=F)
			precocial.matches<-match(precocial.subsample,ecologies$X)
			precocial.variety[which(precocial.variety[,1]==i)[j],2]<-length(which(colSums(ecologies[precocial.matches,3:(dim(ecologies)[2]-1)])>0))
	
			altricial.subsample<-sample(taxa$altricials,i, replace=F)
			altricial.matches<-match(altricial.subsample,ecologies$X)
			altricial.variety[which(altricial.variety[,1]==i)[j],2]<-length(which(colSums(ecologies[altricial.matches,3:(dim(ecologies)[2]-1)])>0))
		}
	}
	return(list(precocial=precocial.variety,altricial=altricial.variety))
}


# The following function will take the estimated ecological diversities by sampling size and fit quantile ribbons  representing
# a 2 sigma confidence interval of ecological richness

make.ribbons<-function(data){
	ribbons<-matrix(NA,1,4)
	colnames(ribbons)<-c('developmental.mode','sample.size','lower','upper')
	ribbons<-data.frame(ribbons)

	for(i in 1:max(data[[1]][,1]) ){
		if(i==1){
			ribbons[1,] <- c('precocial', i, quantile(data$precocial[which(data$precocial[,1]==i),2], probs=c(0.0225, 0.9775) ))
			temporary <- c('altricial',i, quantile(data$altricial[which(data$altricial[,1]==i),2], probs=c(0.0225, 0.9775) ))
			ribbons<-rbind(ribbons,temporary)
		}else{
			temporary <- c('precocial',i, quantile(data$precocial[which(data$precocial[,1]==i),2], probs=c(0.0225, 0.9775) ))
			ribbons<-rbind(ribbons,temporary)

			temporary <- c('altricial',i, quantile(data$altricial[which(data$altricial[,1]==i),2], probs=c(0.0225, 0.9775) ))
			ribbons<-rbind(ribbons,temporary)
		}
	}
	ribbons$sample.size<-as.numeric(ribbons$sample.size)
	ribbons$upper<-as.numeric(ribbons$upper)
	ribbons$lower<-as.numeric(ribbons$lower)
	return(ribbons)
}


# Perform the analyses by applying the functions: 
Bjarnason.flight.rarified <- sample.rarify(taxa=Bjarnason,ecologies=flight.styles,resample.size=1000)
Bjarnason.foot.rarified <- sample.rarify(taxa=Bjarnason,ecologies=footuse.styles,resample.size=1000)
Bj.flight.ribbons<-make.ribbons(Bjarnason.flight.rarified)
Bj.foot.ribbons<-make.ribbons(Bjarnason.foot.rarified)





library(sf)
# This is a package containing functions to facilitate the computation of polygon areas

# This is a function that computes the overlap between polygons representing the richness accumulation curves
# bound by 2-sigma confidence intervals for altricial and precocial populations.
# I am unhappy with the spirit of this approach but I think that it is a bodge that works. 
# We will compare the real overlap coefficients to those expected under scenarios of random ecological category brownian evolution.

compute.overlap <- function(ribbons){
	pre.pol<- st_as_sf(
	rbind(
	setNames(ribbons[which(ribbons[,1]=='precocial'),c(2,4)],c('x','y')),
	apply(setNames(ribbons[which(ribbons[,1]=='precocial'),c(2,3)],c('x','y')),2,rev),
	setNames(ribbons[which(ribbons[,1]=='precocial'),c(2,4)][1,],c('x','y'))
	), coords=c('x','y'))

	alt.pol<- st_as_sf(
	rbind(
	setNames(ribbons[which(ribbons[,1]=='altricial'),c(2,4)],c('x','y')),
	apply(setNames(ribbons[which(ribbons[,1]=='altricial'),c(2,3)],c('x','y')),2,rev),
	setNames(ribbons[which(ribbons[,1]=='altricial'),c(2,4)][1,],c('x','y'))
	), coords=c('x','y'))

	overlap<-st_intersection( sf::st_buffer(st_cast( st_combine(pre.pol$geometry),'POLYGON'), 0), sf::st_buffer(st_cast( st_combine(alt.pol$geometry),'POLYGON'), 0))

	#overlap <- (sf::st_difference(st_cast( st_combine(pre.pol$geometry),'POLYGON'), st_cast( st_combine(alt.pol$geometry),'POLYGON') ))

	fraction.overlap <- sf::st_area(overlap) / ( sf::st_area(overlap)+ (sf::st_area(st_cast( st_combine(pre.pol$geometry),'POLYGON'))- sf::st_area(overlap))+
	(sf::st_area(st_cast( st_combine(alt.pol$geometry),'POLYGON'))- sf::st_area(overlap)) )

	plot(pre.pol,type='l',col='black')
	plot(alt.pol,add=T,col='blue',type='l')
	plot(overlap,add=T,col='green')

	return(fraction.overlap)

}

load('tree.22.10.2022.RData')
# Load the pruned Bjarnason tree

pruned.tree <- keep.tip(pruned.tree, c(Bjarnason$precocials,Bjarnason$altricials)) 
# Prune the tree further to just the taxa of interest

null.flight.overlap<-list()
null.foot.overlap<-list()

# The following loop produces Brownian walks and then ranks and discretises the results 
# to produce fascimiles of ecological trait evolution under a random stochastic process
# We will simulate 1000* possibilities and pass them through the data analytical pipeline
# that was applied to the real data, to determine whether real and simulated datasets
# possess distinct properties (such as whether altricial and precocial groupings determine richness)
# A key component of this approach is that Brownian walk ecology simulations are conducted across
# the same avian phylogeny with the same altricial and precocial groupings as the real data.
# Hence our analysis allows us to parse-out phylogenetic autocorrelation effects that might
# affect the real dataset. 
# Really, 100 simulations would be sufficient. 
# The Brownian walk evolution deciding ecological traits assumes all traits are independent of one another, 
# but the ranking step assumes that the same number of species will develop the trait as is observed in the 
# real population. 
# We therefore simulate a realistic frequency abundance distribution, but not autocorrelation
# among ecological categories. 
# We should think about whether accommodating an ecological autocorrelation structure is necessary but 
# other deficiencies in the method are so large I don't think this degree of finesse is appropriate. 

g<-1
while( g <= 1000){
	brown.flight.styles <- matrix(0,length(flight.styles$X),11) 
	for(i in 1:dim(brown.flight.styles)[2] ){
		uppercut <- length(which(flight.styles[,i+2]==1))
		sim <- fastBM(pruned.tree)
		brown.flight.styles[,i][match(names(sim)[order(sim)[1:uppercut]],flight.styles$X )] <-1 
	}
	brown.flight.styles <- cbind(flight.styles$X,brown.flight.styles) 
	brown.flight.styles <- as.data.frame(brown.flight.styles)
	for(i in 2:12){
		brown.flight.styles[,i] <- as.numeric(brown.flight.styles[,i])
	}
	colnames(brown.flight.styles)[1]<-'X'


	brown.foot.styles <- matrix(0,length(footuse.styles$X),15) 

	for(i in 1:dim(brown.foot.styles)[2] ){
		uppercut <- length(which(footuse.styles[,i+2]==1))
		sim <- fastBM(pruned.tree)
		brown.foot.styles[,i][match(names(sim)[order(sim)[1:uppercut]],footuse.styles$X )] <-1 
	}
	brown.foot.styles <- cbind(footuse.styles$X,brown.foot.styles) 
	brown.foot.styles <- as.data.frame(brown.foot.styles)
	for(i in 2:16){
		brown.foot.styles[,i] <- as.numeric(brown.foot.styles[,i])
	}
	colnames(brown.foot.styles)[1]<-'X'

	Bjarnason.brown.flight.rarified <- sample.rarify(taxa=Bjarnason,ecologies=brown.flight.styles,resample.size=100)
	Bjarnason.brown.foot.rarified <- sample.rarify(taxa=Bjarnason,ecologies=brown.foot.styles,resample.size=100)
	Bj.brown.flight.ribbons<-make.ribbons(Bjarnason.brown.flight.rarified)
	Bj.brown.foot.ribbons<-make.ribbons(Bjarnason.brown.foot.rarified)


	null.flight.overlap[[g]]<- tryCatch({compute.overlap(ribbons=Bj.brown.flight.ribbons)}, error = function(e){},finally={})
	null.foot.overlap[[g]]<- tryCatch({compute.overlap(ribbons=Bj.brown.foot.ribbons)}, error = function(e){},finally={})
	print(g)

	if(is.null(tryCatch({compute.overlap(ribbons=Bj.brown.flight.ribbons)}, error = function(e){},finally={}))==F &
	is.null(tryCatch({compute.overlap(ribbons=Bj.brown.foot.ribbons)}, error = function(e){},finally={})) ==F){
		g<-g+1
	}
}


real.flight.overlap <- compute.overlap(ribbons=Bj.flight.ribbons) 
real.foot.overlap <- compute.overlap(ribbons=Bj.foot.ribbons) 
# What is the overlap between ribbons describing the 2-sigma confidence intervals around ecological richness accumulation
# curves in the altricial and precocial groups? 

length(which(unlist(null.flight.overlap) < real.flight.overlap))/length(unlist(null.flight.overlap))
length(which(unlist(null.foot.overlap) < real.foot.overlap))/length(unlist(null.foot.overlap))
# What proportion of the simulated scenarios involve a greater partition of the precocial and altricial richness curves
# than the real data set? (This is an empirical 1-sided p-value). 


# We will now move on to plotting:

dendr <- dendro_data(as.dendrogram(force.ultrametric(pruned.tree))) 
den.tree <- as.dendrogram(force.ultrametric(pruned.tree))
seg.tree<-dendro_data(den.tree)$segments
# Prepare the phylogeny for phylogram plot production
# The phylogeny differed slightly from an ultrametric topology.

lab.dat<-dendro_data(den.tree)$labels
# Label data

cbp <-c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Palette of colour-blind friendly colours. 

col <- rep(cbp[3],length(lab.dat$label))
col[match(Bjarnason$altricials,lab.dat$label)] <- cbp[2]
# Colour vector of developmental modes 


# Setting colours for the branches to indicate likely precocial and altricial avian lineages: 

branch.col <- rep('Precocial',dim(dendr$segments)[1])
branch.col[which(dendr$segments$x >= 52 & dendr$segments$x <=53 & dendr$segments$xend >= 52 & dendr$segments$xend  <=53 & dendr$segments$y  <=132)] <- 'Altricial'
branch.col[which(dendr$segments$x >= 47 & dendr$segments$x <=50 & dendr$segments$xend >= 47 & dendr$segments$xend  <=50 & dendr$segments$y  <=100)] <- 'Altricial'
branch.col[which(dendr$segments$x == 29  & dendr$segments$xend == 29  & dendr$segments$y  <=150)] <- 'Altricial'
branch.col[which(dendr$segments$x >= 20 & dendr$segments$x <=27 & dendr$segments$xend >= 20 & dendr$segments$xend  <=27 & dendr$segments$y  <=124)] <- 'Altricial'
branch.col[which(dendr$segments$x >= 1 & dendr$segments$x <=19 & dendr$segments$xend >= 1 & dendr$segments$xend  <=19 & dendr$segments$y  <=125)] <- 'Altricial'

branch.col[which(dendr$segments$x >= 54 & dendr$segments$x <=68 & dendr$segments$xend >= 54 & dendr$segments$xend  <=68 & dendr$segments$y  <=123)] <- 'Precocial'
branch.col[which(dendr$segments$x == 51 & dendr$segments$x ==51 & dendr$segments$xend == 51 & dendr$segments$xend  ==51 & dendr$segments$y  <=127)] <- 'Precocial'
branch.col[which(dendr$segments$x >= 30 & dendr$segments$x <=46 & dendr$segments$xend >= 30 & dendr$segments$xend  <=46 & dendr$segments$y  <=125)] <- 'Precocial'
branch.col[which(dendr$segments$x == 28 & dendr$segments$x ==28 & dendr$segments$xend == 28 & dendr$segments$xend  ==28 & dendr$segments$y  <=130)] <- 'Precocial'
df.phy <- segment(dendr)
df.phy <- cbind(df.phy,branch.col)

# Information for labelling segments of the phylogeny by their major subclade/grade:

Passeriformes <- c(1,5)
Path_Passeriformes <- data.frame(x=Passeriformes, y=c(-2, -2))

Psittaciformes <- c(6,9)
Path_Psittaciformes <- data.frame(x=Psittaciformes, y=c(-2, -2))

Coraciimorphae <- c(10,17)
Path_Coraciimorphae <- data.frame(x=Coraciimorphae, y=c(-2, -2))

Strigiformes <- c(18,19)
Path_Strigiformes <- data.frame(x=Strigiformes, y=c(-2, -2))

Aequorlitornithes <- c(20,46)
Path_Aequorlitornithes <- data.frame(x=Aequorlitornithes, y=c(-2, -2))

Columbaves <- c(47,51)
Path_Columbaves <- data.frame(x=Columbaves, y=c(-2, -2))

Apodiformes <- c(52,53)
Path_Apodiformes <- data.frame(x=Apodiformes, y=c(-2, -2))

Galloanserae <- c(54,66)
Path_Galloanserae <- data.frame(x=Galloanserae, y=c(-2, -2))

Galloanserae <- c(54,66)
Path_Galloanserae <- data.frame(x=Galloanserae, y=c(-2, -2))

Palaeognathae <- c(67,68)
Path_Palaeognathae <- data.frame(x=Palaeognathae, y=c(-2, -2))


library(geomtextpath) # version 0.1.1
# Package to bend text around specified paths. 

binary.flight.scores <- flight.styles[,3:13] 
binary.flight.scores<-as.matrix(binary.flight.scores)
# Prepare the data as a matrix to make it easy to index. 
binary.flight.scores[which(binary.flight.scores=='?' | binary.flight.scores=='' )] <- 0 
binary.flight.scores[grep(binary.flight.scores,pattern='\\?')] <- 0 
binary.flight.scores <- apply(binary.flight.scores,2,FUN=as.numeric)
# I am making sure the data class is numeric
rownames(binary.flight.scores)<-flight.styles$X
# Set the row names of the flight style matrix to our taxa 

binary.foot.scores <- footuse.styles[,3:17] 
binary.foot.scores<-as.matrix(binary.foot.scores)
# Prepare the data as a matrix to make it easy to index. 
binary.foot.scores[which(binary.foot.scores=='?' | binary.foot.scores=='' )] <- 0 
binary.foot.scores[grep(binary.foot.scores,pattern='\\?')] <- 0 
binary.foot.scores <- apply(binary.foot.scores,2,FUN=as.numeric)
# I am making sure the data class is numeric
rownames(binary.foot.scores)<-footuse.styles$X
# Set the row names of the foot style matrix to our taxa 

binary.flight.scores<-binary.flight.scores[pruned.tree$tip,]
df <- data.frame(binary.flight.scores)
# coerce matric to dataframe in preparation for plotting 
# coerce matric to dataframe in preparation for plotting 
generic <- rep('0',dim(df)[1])

generic[which(rowSums(df)==0)]<-'1'
generic<-as.numeric(generic)
df<-cbind(df,generic)
# Species which were not scores as having any particular flight speciality. 

colnames(df)[12] <- 'no identified specialty'
# Label the column with no particular flight speciality. 

library(reshape2)
# Package for re-arranging data into different data table formats

df2<-cbind(rownames(df),df)
colnames(df2)[1]<-'taxon'
df2<-melt(df2)
df2<-df2[-which(df2$value=='0'),]
df2$variable <- factor(df2$variable, levels= c(names(table(df2$variable))[order(table(df2$variable))]))
# Data frame of all considered bird species and all of their flight styles. 

lab.order<-dendro_data(den.tree)$labels
df2$taxon <- factor(df2$taxon, levels= lab.order$label)
df2$variable <- factor(df2$variable, levels= c(names(table(df2$variable))[order(table(df2$variable))]))
# Render the data as factor variables. 

scale <- 0.1
size<-4/.pt
yspace<- 0
# Plot parameters 

 Subplot.title.font.size <- 30/.pt
 Subplot.axis.font.size <- 16/.pt
 Legend.axis.font.size <- 16/.pt
 Legend.title.font.size <- 16/.pt
 Legend.tick.size <- 1
 Subplot.linewidth <- 1/2
 Subplot.tick.size <-1
 legend.key.width <- 1
 legend.key.height <- 0.5

# Produce a phylogram of flight styles across the study birds, 
# flight ecological categories are represented as stacked bar charts, 
# and lineages are coloured by their likely developmental strategy affiliation
# (in reality a variety of hypotheses might describe variation in developmental strategy across this tree, 
# but the best supported versions will resemble this)

phy.plot.bjar.flight <- 
ggplot()+
	geom_bar(aes( x=taxon, fill=variable), data= df2,position= position_stack(reverse=F, vjust=0))+
	scale_fill_manual(values= c("#D55E00","#F0E442", "#FFFF00", "#0072B2", "#CC0000","#000000", "#FFFFFF", "#009E73", "#56B4E9", 
                     "#CC79A7", "#FF0000", "#CC00CC", "#0000FF"), 
	labels=c('Flap-bounding', 'Hovering', 'Terrestrial soaring', 'Sallying','Wing-propelled diving','Aerial pursuit hunting',
	'no identified specialty','Ocean soaring','Escape burst','Flap-gliding','Migrant','Cluttered') )+
	geom_segment(
	data=df.phy,
  	aes(x = x, y = (-y *scale)-1, xend = xend, yend = (-yend *scale)-1, col=branch.col), 
      lineend = "square", lwd=Subplot.linewidth)+
	scale_colour_manual(values=c(cbp[2],cbp[3]))+
	guides(colour='none')+
	labs(fill='Flight-style')+
	scale_y_continuous(expand = c(0.2, 0), labels = function(y) y + 1) + 
	geom_rect(aes(xmin=0,xmax=max(df.phy$x)+0.1,ymin=-1,ymax=0), fill='white',alpha=0.5)+
	geom_rect(aes(xmin=0,xmax=max(df.phy$x)+0.1,ymin=-2,ymax=0), fill='white',alpha=0.5)+
	geom_rect(aes(xmin=0,xmax=max(df.phy$x)+0.1,ymin=-5,ymax=0), fill='white',alpha=0.5)+
	geom_rect(aes(xmin=0,xmax=max(df.phy$x)+0.1,ymin=-10,ymax=0), fill='white',alpha=0.5)+
	geom_textpath(data=Path_Passeriformes,size  = size, label = 'Passeriformes', aes(x=x,y=y-yspace), linecolour=NA)+
	geom_segment(data= Path_Passeriformes,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='darkgrey')+ 
	geom_textpath(data=Path_Psittaciformes,size  = size, label = 'Psittaciformes', aes(x=c(6.75,9.75),y=y-yspace), linecolour=NA)+
	geom_segment(data= Path_Psittaciformes,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='black')+ 
	geom_textpath(data=Path_Coraciimorphae,size  = size, label = 'Coraciimorphae', aes(x=x,y=y-yspace), linecolour=NA)+
	geom_segment(data= Path_Coraciimorphae,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='darkgrey')+ 
	geom_segment(data= Path_Strigiformes,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='black')+ 
	geom_textpath(data=Path_Strigiformes,size  = size, label = 'Strigiformes', aes(x=c(16,21),y=y-yspace), linecolour=NA)+
	geom_textpath(data=Path_Aequorlitornithes,size  = size, label = 'Aequorlitornithes', aes(x=x,y=y-yspace), linecolour=NA)+
	geom_segment(data= Path_Aequorlitornithes,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='darkgrey')+ 
	geom_textpath(data=Path_Columbaves,size  = size, label = 'Columbaves', aes(x=x-1,y=y-yspace), linecolour=NA)+
	geom_segment(data= Path_Columbaves,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='black')+ 
	geom_textpath(data=Path_Apodiformes,size  = size, label = 'Apodiformes', aes(x=c(50,55),y=y-yspace), linecolour=NA)+
	geom_segment(data= Path_Apodiformes,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='darkgrey')+ 
	geom_textpath(data=Path_Galloanserae,size  = size, label = 'Galloanserae', aes(x=x,y=y-yspace), linecolour=NA)+
	geom_segment(data= Path_Galloanserae,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='black')+ 
	geom_textpath(data=Path_Palaeognathae,size  = size, label = 'Palaeognathae', aes(x=x+c(-2.5,-.5),y=y-yspace), linecolour=NA)+
	geom_segment(data= Path_Palaeognathae,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='darkgrey')+ 
	#geom_text(data= lab.dat, aes(x=x, y=y, label=gsub('_.*','',label)), angle=(createAngleHJustCols(lab.dat)[['angle']])+(seq(-5/2,5/2,length.out=68)), hjust=(createAngleHJustCols(lab.dat)[['hjust']]),size = 12 / .pt, col='black' )+
	coord_polar()+
	#lims(y=c(20,-10))+
	expand_limits(x=0)+
	theme(legend.position=c(1.00,0.5),
	#legend.background=element_blank(),
	#legend.title=element_blank(),
	legend.text=element_text(size=Legend.axis.font.size),
	legend.title=element_text(size=Legend.axis.font.size),
	legend.key.height=unit(0.1,'cm'),
	legend.key.width=unit(0.1,'cm'),
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	axis.ticks.x=element_blank(),
	plot.margin=margin(t =  -25, # Top margin
                             r = -40, # Right margin
                             b = -45, # Bottom margin
                            l = -110) # Left margin
      )


# We will now repeat the same approach for foot-use scores across the study birds. 

binary.foot.scores<-binary.foot.scores[pruned.tree$tip,]
df <- data.frame(binary.foot.scores)
# coerce matrix to dataframe in preparation for plotting 

df2<-cbind(rownames(df),df)
colnames(df2)[1]<-'taxon'
df2<-melt(df2)
df2<-df2[-which(df2$value=='0'),]
df2$variable <- factor(df2$variable, levels= c(names(table(df2$variable))[order(table(df2$variable))]))

lab.order<-dendro_data(den.tree)$labels
df2$taxon <- factor(df2$taxon, levels= lab.order$label)
df2$variable <- factor(df2$variable, levels= c(names(table(df2$variable))[order(table(df2$variable))]))
# Coerce data to factor class

# Produce a phylogram of foot-use variety across the precocial and altricial birds
phy.plot.bjar.foot<-
ggplot()+
	geom_bar(aes( x=taxon, fill=variable), data= df2,position= position_stack(reverse=F, vjust=0))+
	scale_fill_manual(values= c("#F0E442","#009E73", "#0072B2", "#CC00CC", "#D55E00","#000000", "#CC79A7", "#56B4E9", "#FF0000", 
                     "#0000FF", "#FFFF00", "#CC0000", "#E69F00",'darkgrey','grey'), 
	labels=c('Scansorial','Raptorial','Digit manipulation','Hopping','Forceful grip','Holding down','Hanging clinging',
	'Submerged foot propelled swimming','Ground running specialist','Surface swimming','Obligate walking','Swimming ability',
	'Perching','Documented walking ability','Facultative walking') )+
	geom_segment(
	data=df.phy,
  	aes(x = x, y = (-y *scale)-1, xend = xend, yend = (-yend *scale)-1, col=branch.col), 
      lineend = "square", lwd= Subplot.linewidth)+
	scale_colour_manual(values=c(cbp[2],cbp[3]))+
	guides(colour='none')+
	labs(fill='Foot-use')+
	scale_y_continuous(expand = c(0.2, 0), labels = function(y) y + 1) + 
	geom_rect(aes(xmin=0,xmax=max(df.phy$x)+0.1,ymin=-1,ymax=0), fill='white',alpha=0.5)+
	geom_rect(aes(xmin=0,xmax=max(df.phy$x)+0.1,ymin=-2,ymax=0), fill='white',alpha=0.5)+
	geom_rect(aes(xmin=0,xmax=max(df.phy$x)+0.1,ymin=-5,ymax=0), fill='white',alpha=0.5)+
	geom_rect(aes(xmin=0,xmax=max(df.phy$x)+0.1,ymin=-10,ymax=0), fill='white',alpha=0.5)+
	geom_textpath(data=Path_Passeriformes,size  = size, label = 'Passeriformes', aes(x=x,y=y-yspace), linecolour=NA)+
	geom_segment(data= Path_Passeriformes,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='darkgrey')+ 
	geom_textpath(data=Path_Psittaciformes,size  = size, label = 'Psittaciformes', aes(x=c(6.75,9.75),y=y-yspace), linecolour=NA)+
	geom_segment(data= Path_Psittaciformes,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='black')+ 
	geom_textpath(data=Path_Coraciimorphae,size  = size, label = 'Coraciimorphae', aes(x=x,y=y-yspace), linecolour=NA)+
	geom_segment(data= Path_Coraciimorphae,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='darkgrey')+ 
	geom_segment(data= Path_Strigiformes,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='black')+ 
	geom_textpath(data=Path_Strigiformes,size  = size, label = 'Strigiformes', aes(x=c(16,21),y=y-yspace), linecolour=NA)+
	geom_textpath(data=Path_Aequorlitornithes,size  = size, label = 'Aequorlitornithes', aes(x=x,y=y-yspace), linecolour=NA)+
	geom_segment(data= Path_Aequorlitornithes,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='darkgrey')+ 
	geom_textpath(data=Path_Columbaves,size  = size, label = 'Columbaves', aes(x=x-1,y=y-yspace), linecolour=NA)+
	geom_segment(data= Path_Columbaves,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='black')+ 
	geom_textpath(data=Path_Apodiformes,size  = size, label = 'Apodiformes', aes(x=c(50,55),y=y-yspace), linecolour=NA)+
	geom_segment(data= Path_Apodiformes,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='darkgrey')+ 
	geom_textpath(data=Path_Galloanserae,size  = size, label = 'Galloanserae', aes(x=x,y=y-yspace), linecolour=NA)+
	geom_segment(data= Path_Galloanserae,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='black')+ 
	geom_textpath(data=Path_Palaeognathae,size  = size, label = 'Palaeognathae', aes(x=x+c(-2.5,-.5),y=y-yspace), linecolour=NA)+
	geom_segment(data= Path_Palaeognathae,aes(x=x[1]-.5, xend=x[2]+.5, y=-1, yend=-1), linewidth=Subplot.linewidth, colour='darkgrey')+ 
	#geom_text(data= lab.dat, aes(x=x, y=y, label=gsub('_.*','',label)), angle=(createAngleHJustCols(lab.dat)[['angle']])+(seq(-5/2,5/2,length.out=68)), hjust=(createAngleHJustCols(lab.dat)[['hjust']]),size = 12 / .pt, col='black' )+
	coord_polar()+
	#lims(y=c(20,-10))+
	expand_limits(x=0)+
	theme(legend.position=c(1.05,0.5),
	#legend.background=element_blank(),
	#legend.title=element_blank(),
	legend.text=element_text(size=Legend.axis.font.size),
	legend.title=element_text(size=Legend.axis.font.size),
	legend.key.height=unit(0.1,'cm'),
	legend.key.width=unit(0.1,'cm'),
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	axis.ticks.x=element_blank(),
	plot.margin=margin(t =  -48, # Top margin
                             r = -40,# Right margin
                             b = -40,# Bottom margin
                            l = -120) # Left margin
      )


# Now it is time to plot the richness accumulation curves 

# Flight-style
bj.flight<-
ggplot(data=Bj.flight.ribbons)+
geom_ribbon(aes(ymin=lower,ymax=upper,x=sample.size, colour=developmental.mode),lwd= Subplot.linewidth ,alpha=.5, fill=NA)+
scale_colour_manual(values=c(cbp[2],cbp[3]))+
labs(x='Sample size',y='Number of activities', fill='developmental mode', colour='developmental mode')+
geom_ribbon(aes(fill=developmental.mode,ymin=lower,ymax=upper,x=sample.size, colour=developmental.mode),alpha=.5,lwd= Subplot.linewidth)+
scale_fill_manual(values=c(cbp[2],cbp[3]))+
labs(x='Sample size',y='Number of activities')+
ggtitle('Flight-style variety')+
theme(
legend.position='bottom')+
  theme(axis.line.x.bottom=element_line(size = Subplot.linewidth, colour = "black"),
	axis.line.y.left=element_line(size = Subplot.linewidth, colour = "black"),
      axis.text.x=element_text(size= Subplot.axis.font.size,colour='black'),
      axis.text.y=element_text(size= Subplot.axis.font.size,colour='black'),
	#axis.title.x.bottom=element_text(size=16,colour='black'),
	#axis.title.y.left=element_text(size=16,colour='black', fill='developmental mode', colour='developmental mode'),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
      axis.ticks=element_line(size=Subplot.tick.size),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.title=element_text(size=Legend.axis.font.size,colour='black'),
	legend.text=element_text(size=Legend.axis.font.size,colour='black'),
	plot.title=element_text(size=Subplot.title.font.size,colour='black',hjust=0.5))


# Foot-use
bj.foot<-
ggplot(data=Bj.foot.ribbons)+
geom_ribbon(aes(ymin=lower,ymax=upper,x=sample.size, colour=developmental.mode),lwd= Subplot.linewidth ,alpha=.5, fill=NA)+
scale_colour_manual(values=c(cbp[2],cbp[3]))+
labs(x='Sample size',y='Number of activities', fill='developmental mode', colour='developmental mode')+
geom_ribbon(aes(fill=developmental.mode,ymin=lower,ymax=upper,x=sample.size, colour=developmental.mode),alpha=.5,lwd= Subplot.linewidth)+
scale_fill_manual(values=c(cbp[2],cbp[3]))+
labs(x='Sample size',y='Number of activities')+
ggtitle('Foot-use variety')+
theme(
legend.position='bottom')+
  theme(axis.line.x.bottom=element_line(size = Subplot.linewidth, colour = "black"),
	axis.line.y.left=element_line(size = Subplot.linewidth, colour = "black"),
      axis.text.x=element_text(size= Subplot.axis.font.size,colour='black'),
      axis.text.y=element_text(size= Subplot.axis.font.size,colour='black'),
	#axis.title.x.bottom=element_text(size=16,colour='black'),
	#axis.title.y.left=element_text(size=16,colour='black', fill='developmental mode', colour='developmental mode'),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
      axis.ticks=element_line(size=Subplot.tick.size),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	legend.title=element_text(size=Legend.axis.font.size,colour='black'),
	legend.text=element_text(size=Legend.axis.font.size,colour='black'),
	plot.title=element_text(size=Subplot.title.font.size,colour='black',hjust=0.5))



library(ggpubr)
# For arranging plots
library(grid)
# Package for arranging plots. 

top <- ggarrange(bj.flight,bj.foot,ncol=1,labels=c('b','d'), vjust=1, font.label = list(size = 20, color = "black", face = "bold", family = NULL),
common.legend=T, legend='top')
library(grid)
top <- annotate_figure(top, left = textGrob("Richness", rot = 90, vjust = 1, gp = gpar(cex = 1.3,fontsize=30/.pt)),
                    bottom = textGrob("Sample size", gp = gpar(cex = 1.3,fontsize=30/.pt), vjust=.1))

phys <- ggarrange( phy.plot.bjar.flight, phy.plot.bjar.foot, ncol=1,labels=c('a','c'), vjust=c(2.2,1), font.label = list(size = 20, color = "black", face = "bold", family = NULL ))
ggarrange(phys,top,ncol=2,widths=c(1,.8))
# Compile and label subplots 


#setwd('C:/Users/Lab/Documents/AOrkney/Precocial_study')
setwd('Z:/Andrew_backup/Precocial_paper')
ggsave(filename='FigII_10_30_2024.pdf',height=13,width=17.8,unit='cm', dpi=300) 
# Save plot


