# This is a script that produces phenograms/ 'bird nest plots' of divergent wing and leg size evolution across Bjarnason's
# birds, under a two-rate Brownian motion model depending on developmental strategy.
# The phenogram branches are decorated with colours indicating likely developmental strategy, terminal branches that
# lead to altricial species practicing novel flight-styles or foot-uses are coloured black.
# Histograms of flight-styles and foot-uses that only occur in altricial birds will decorate the time=0 horizon. 
# The time=0 horizon will also be decorated with silhouettes illustrating a range of phenotypes of altricial birds across
# trait space. 
#
# This script is waiting for new silhouettes that will be contributed by Isaiah and Will 

library( ape ) # 5.7-1
library(phytools) # 1.5-1
library( geomorph ) # 4.0.5
library(mvMORPH) # 1.1.8

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
# Render developmental strategy data consistent 

tips <- names(masses)[which(developmental_mode == 'Precocial' | developmental_mode == 'Altricial')]
# Define taxa with well-defined precocial and altricial developmental strategies 

short.tree <- keep.tip(pruned.tree,tips)
short.mode <- developmental_mode[short.tree$tip]
# Prune data to taxa of interest

tree<-make.simmap(short.tree,short.mode , model="ER", nsim=1, pi= c(0,1) )
# Change order of pi to assume a different root state
col<-c("blue","orange"); names(col)<-c("Precocial","Altricial")
plotSimmap(tree,col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)
# Construct a single possible evolutionary history of developmental strategy, treated as a stochastic character
# Precocial development is assumed to be ancestral

leg <- (GPA.Csize$tibiotarsus+GPA.Csize$femur+GPA.Csize$tarsometatarsus)
wing <- (GPA.Csize$humerus+GPA.Csize$radius+GPA.Csize$carpometacarpus)
total <- wing+leg
trait <- ((wing-leg)/total)
# Define the trait as the difference between the wing and leg, as compared to their total size
# This is a measure of limb skeleton 'evenness'

model <- mvBM(tree=tree, data= trait[tree$tip] , model='BMM')
# Fit a single 2-rate Brownian motion model of trait evolution depending on developmental strategy across the tree

#setwd('C:/Users/Lab/Documents/AOrkney/bird_scripts_and_data_15_02_2023')
Foot.use <- read.csv('standardised_foot_scores_11_07_2023.csv')
Flight.style <- read.csv('flight_masses_22_10_2022_plus_A.csv')
# Load ecological classifications for Bjarnason's birds

# Identify foot-uses that are specific to altricial birds in Bjarnason's dataset
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

# Perform the same operation seeking altricial-specific flight-styles
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

library( ggplot2 ) # 3.4.1
library( ggdendro ) # 0.1.23
library( dendextend ) # 1.17.1
library(zoo) # 1.8-12
library(dplyr) # 1.1.1

# The following function is deployed to construct a phenogram of ancestral state reconstructions, 
# across the bird family tree, under a specific model. The mvMORPH 'estim' function is integral to this process.

make.dendr <- function(tree, trait, mod ){
	dendr <- dendro_data(as.dendrogram(force.ultrametric(tree)))
	lab.dat <- dendr$labels

	fit<-estim(tree, data= trait, object=mod, asr=TRUE)
	Anc <- rep(NA, length(dendr$segments[,1])) 
	Anc [which(dendr$segments$x==dendr$segments$xend)] <- fit$estim[match(tree$edge[,1],rownames(fit$estim))] 
	Anc [which(dendr$segments$x==dendr$segments$xend)-1] <- fit$estim[match(tree$edge[,1],rownames(fit$estim))]

	dendr.mod<-dendr$segments/2
	dendr.mod$x[which(dendr$segments$x==dendr$segments$xend)-1] <- fit$estim[match(tree$edge[,1],rownames(fit$estim))] # This defines the nodal positions of the nodes
	dendr.mod$xend[which(dendr$segments$x==dendr$segments$xend)-1] <- dendr.mod$x[which(dendr$segments$x==dendr$segments$xend)-1]# This collapses nodes by defining their ends as the same as their starts
	dendr.mod$x[which(dendr$segments$x==dendr$segments$xend)] <-   fit$estim[match(tree$edge[,1],rownames(fit$estim))] # This defines the nodal positions of the start of each branch
	dendr.mod$xend[which(dendr$segments$yend==0)] <- trait[match(lab.dat$label,names(trait)) ] # This defines the final position of each terminal branch
	
	for(i in 1:dim(dendr.mod)[1]){
		if(dendr$segments[i,]$yend !=0 ){ # We do not want to move distal tips of branches that hit zero
			if( is.na(match(dendr.mod[i,]$xend,dendr.mod$x))==T  ){ # We check that the branch has no descendants 
				choice <- which(dendr$segments$x==dendr$segments[i,]$xend) # We find a potential site it should attach to 
				target <- choice[which(dendr$segments[choice,]$y-dendr$segments[choice,]$yend == 0)][1] # We have found a node to target
				dendr.mod$xend[i]<- dendr.mod$x[target]
			
			}
		}	
	}

	Rate <- abs((dendr.mod$xend-dendr.mod$x))/(dendr.mod$y-dendr.mod$yend)

	return( cbind(dendr.mod, Anc,Rate ) )

}

dendr.trait <- make.dendr(tree=tree, trait=trait[short.tree$tip], mod=model)
# Produce an ancestral state reconstruction of trait evolution in the context of a two-rate model that depends on developmental mode

dendr <- dendro_data(as.dendrogram(force.ultrametric(short.tree)))
lab.dat <- dendr$labels
# Compute the phylogeny structure format and labels necessary for plotting in ggplot2


cbp <-c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Define a palette of colours that are friendly to most human visual differences

col <- rep(cbp[3],length(lab.dat$label))
col[match(names(developmental_mode[which(developmental_mode=='Altricial')]),lab.dat$label)] <- cbp[2]
# Define a dichotomous vector of colours to denote altricial and precocial development across the tips of the family tree

	df.phy <- segment(dendr)
	branch.col <- rep('Precocial',dim(dendr$segments)[1])
	branch.col[which(dendr$segments$x >= 1 & dendr$segments$x <=60 & dendr$segments$xend >= 1 & dendr$segments$xend  <=60 & dendr$segments$y  <=125)] <- 'Altricial'
	branch.col[which(dendr$segments$x >= 61 & dendr$segments$x <=70 & dendr$segments$xend >= 61 & dendr$segments$xend  <=70 & dendr$segments$y  <=125)] <- 'Altricial'
	branch.col[which(dendr$segments$x ==72 & dendr$segments$xend ==72 & dendr$segments$y  <=135)] <- 'Altricial'
	branch.col[which(dendr$segments$x >= 92 & dendr$segments$x <=96 & dendr$segments$xend >= 92 & dendr$segments$xend  <=96 & dendr$segments$y  <=70)] <- 'Altricial'
	branch.col[which(dendr$segments$x >= 92 & dendr$segments$x <=96 & dendr$segments$xend >= 92 & dendr$segments$xend  <=96 & dendr$segments$y  <=70)] <- 'Altricial'
	branch.col[which(dendr$segments$x ==98 & dendr$segments$xend ==98 & dendr$segments$y  <=145)] <- 'Altricial'
	branch.col[which(dendr$segments$x >= 100 & dendr$segments$x <=104 & dendr$segments$xend >= 100 & dendr$segments$xend  <=104 & dendr$segments$y  <=130)] <- 'Altricial'
	# Colour the individual branches of the family tree with likely developmental strategies; this is basically just a cartoon

createAngleHJustCols <- function(labeldf) {        
    nn <- length(labeldf$y)
    halfn <- floor(nn/2)
    firsthalf <- rev(90 + seq(0,360, length.out = nn))
    secondhalf <- rev(-90 + seq(0,360, length.out = nn))
    angle <- numeric(nn)
    angle[1:halfn] <- firsthalf[1:halfn]
    angle[(halfn+1):nn] <- secondhalf[(halfn+1):nn]

    hjust <- numeric(nn)
    hjust[1:halfn] <- 0
    hjust[(halfn+1):nn] <- 1

    return(list(angle = angle, hjust = hjust))
}
# This is a function that I borrowed from Stackexchange: https://stackoverflow.com/questions/38034663/rotate-labels-for-ggplot-dendrogram
# The function rotates labels on fan phylogenies. 

scale <- 0.1
size<-10/.pt
yspace<-.1
# These are plotting parameters which will need to be changed before the final compilation of a journal-ready figure

# Produce a test plot of the bird family tree with branches coloured by likely developmental strategy
phy.plot <- 
ggplot()+
	geom_segment(
	data=df.phy,
  	aes(x = x, y = y *scale, xend = xend, yend = yend *scale, col=branch.col), 
      lineend = "square", lwd=3/2)+
	scale_colour_manual(values=c(cbp[2],cbp[3]))+
	scale_y_reverse(expand = c(0.2, 0)) + 
	geom_text(data= lab.dat, aes(x=x, y=y, label=gsub('_.*','',label)), angle=(createAngleHJustCols(lab.dat)[['angle']]), hjust=(createAngleHJustCols(lab.dat)[['hjust']]),size = 6 / .pt, col=col )+
	coord_polar()+
	lims(y=c(20,-10))+
	expand_limits(x=0)+
	theme(legend.position=c(0.5,0.5),
	legend.background=element_blank(),
	legend.title=element_blank(),
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
	plot.margin=margin(t = -80,  # Top margin
                             r = -90,  # Right margin
                             b = -90,  # Bottom margin
                             l = -80) # Left margin
      )


dendr.bird.score <- cbind(dendr.trait,branch.col)
# Bind the colour vector to the dendrogram trait object


	dendr.bird.score <- cbind(dendr.bird.score,rep('NA',dim(dendr.trait)[1]))
	colnames(dendr.bird.score)[8]<-'nov'
	vector <- which(dendr$segments$yend==0)
	dendr.bird.score$nov[vector[which(lab.dat$label %in% Novel_flight==T)]]<-'Novel'
	# Identify terminal branches of the bird family tree that lead to birds practicing novel flight-styles
	dendr.novel<-dendr.bird.score
	dendr.novel$branch.col[which(dendr.bird.score$nov=='Novel')]<-'Novel Altricial'

library(reshape2) # 1.4.4 
# For data formatting and management

df<-data.frame(Flight.style[c(3,6,11,13)])
rownames(df)<-Flight.style[,1]
df2<-cbind(rownames(df),df,dendr.bird.score$xend[which(dendr.bird.score$yend==0)][match(rownames(df), dendr$labels$label)] )
df2<-df2[complete.cases(df2),]
colnames(df2)[1]<-'taxon'
colnames(df2)[dim(df2)[2]]<-'ypos'
df2<-melt(df2, id=c('taxon','ypos'))
df2<-df2[-which(df2$value=='0'),]
df2$variable <- factor(df2$variable, levels= c(names(table(df2$variable))[order(table(df2$variable))]))
lab.order<-dendr$labels$label
df2$taxon <- factor(df2$taxon, levels= lab.order)
df2<-rbind(df2,df2)
# Form a data frame to carry the information needed for plotting a histogram distribution of flight-style categorizations
# across trait-space 


df.sub <- df2[,1:2]
df.sub <-df.sub[!duplicated(df.sub), ]
df.sub <- df.sub[order(df.sub$ypos),]
df.sub <- df.sub[c(1,3,4,7, 11,12,14,17, 22, 24, 32,34,38,39,42,43,44,45,46,47,49,50,51),]
# Eliminate duplications

library(png) # 0.1-8
library(grid) # base

#setwd('C:/Users/Lab/Documents/AOrkney/Precocial_study/bird silhouettes')
# Set path to load silhouettes 

#Upupa_epops_silhouette <- readPNG('Upupa_epops_silhouette.png')
#Upupa_epops_silhouette[,,2][which(Upupa_epops_silhouette[,,2]==1)]<-NA
#Upupa_epops_silhouette[,,2][which(Upupa_epops_silhouette[,,2]==0)]<-1
#Upupa_epops_silhouette[,,2][which(Upupa_epops_silhouette[,,2]!=1)]<-NA

#Upupa.raster <- as.raster(Upupa_epops_silhouette[,,2], interpolate=F)
#Upupa.raster[Upupa.raster == "#FFFFFF"] <- "#CCCCCC"
#Upupa.g = rasterGrob(Upupa.raster,interpolate=T)

# Define a function to prepare silhouette images for plotting 
prep.image <- function(string){
	string <- readPNG(string)
	string[,,2][which(string[,,2]==1)]<-NA
	#string[,,2][which(string[,,2]==0)]<-1
	string[,,2][which(string[,,2]<0.1)]<-1
	string[,,2][which(string[,,2]!=1)]<-NA
	string.raster <- as.raster(string[,,2],interpolate=F)
	string.raster[string.raster=="#FFFFFF"] <- "#CCCCCC"
	string.g <- rasterGrob(string.raster,interpolate=T)
	return(string.g)
}

Upupa.g <- prep.image('Upupa_epops_silhouette.png')
Topaza.g <- prep.image('Topaza_pyra_silhouette.png')
Pelecanus.g <- prep.image('Pelecanus_occidentalis_silhouette.png')
Pelagodroma.g <- prep.image('Pelagodroma_marina_silhouette.png')
Nyctibius.g <- prep.image('Nyctibius_griseus_silhouette.png')
Leucocarbo.g <- prep.image('Leucocarbo_atriceps_silhouette.png')
Jynx.g <- prep.image('Jynx_torquilla_silhouette.png')
Fregata.g <- prep.image('Fregata_silhouette.png')
Columba.g <- prep.image('Columba_palumbus_silhouette.png')
# Load and prepare images




# Produce decorated phenogoram plot: 
dev.tree.flight <- 
ggplot()+
	annotation_custom(grob=Upupa.g, xmin=0, xmax=20, ymin=-0.01214379-0.05, ymax=-0.01214379+0.05)+
	annotation_custom(grob=Topaza.g, xmin=0, xmax=20, ymin=-0.28959137-0.05, ymax=-0.28959137+0.05)+
	annotation_custom(grob=Pelecanus.g, xmin=2, xmax=22, ymin=0.32817671-0.05, ymax=0.32817671+0.05)+
	annotation_custom(grob=Pelagodroma.g, xmin=4, xmax=24, ymin=-0.39730291-0.05, ymax=-0.39730291+0.05)+
	annotation_custom(grob=Nyctibius.g, xmin=-6, xmax=16, ymin=0.26653139-0.02, ymax=0.26653139+0.04)+
	annotation_custom(grob=Leucocarbo.g, xmin=-5, xmax=20, ymin=0.10746715-0.01, ymax=0.10746715+0.1)+
	annotation_custom(grob=Jynx.g, xmin=0, xmax=20, ymin=-0.17065767-0.03, ymax=-0.17065767+0.03)+
	annotation_custom(grob=Fregata.g, xmin=0, xmax=20, ymin=0.50721715-0.08, ymax=0.50721715+0.05)+
	annotation_custom(grob=Columba.g, xmin=0, xmax=20, ymin=-0.12040668-0.02, ymax=-0.12040668+0.06)+

	geom_histogram(aes( y=ypos, fill=variable), data= df2, binwidth=0.02, alpha=.5)+
	xlim(-72.907238,19)+
	coord_flip()+
	scale_fill_manual(values= c( "#009E73", 
                      "#F0E442", "#0072B2", "#D55E00", "#CC79A7", 'darkgrey', 'lightgrey'),
	labels=c('Ocean soaring','Hovering','Sallying','Flap-bounding') )+
	geom_segment(data=dendr.bird.score[which(dendr.bird.score$branch.col=='Altricial'),],
  	aes(x = -y, y = x, xend = -yend, yend =xend ),col='black', lwd=3)+
	geom_segment(data=dendr.bird.score[which(dendr.bird.score$branch.col=='Altricial'),],
  	aes(x = -y, y = x, xend = -yend, yend =xend ),col=cbp[2], lwd=3/2)+

	geom_segment(data=dendr.bird.score[which(dendr.bird.score$branch.col=='Precocial'),],
  	aes(x = -y, y = x, xend = -yend, yend =xend ),col='black', lwd=3)+
	geom_segment(data=dendr.bird.score[which(dendr.bird.score$branch.col=='Precocial'),],
  	aes(x = -y, y = x, xend = -yend, yend =xend ),col=cbp[3], lwd=3/2)+
	#labs(col= expression(italic(sigma)), x=expression('Mya'), y=expression('(wing-leg)/(leg+wing)'))+

	geom_segment(data=dendr.bird.score[which(dendr.bird.score$nov=='Novel'),],
  	aes(x = -y, y = x, xend = -yend, yend =xend ),col='black', lwd=3/2)+ 

	geom_segment(aes(x=-40,y=0.4,xend= -40,yend=0.5),col='black',lwd=3)+
	geom_segment(aes(x=-50,y=0.4,xend= -50,yend=0.5),col='black',lwd=3)+
	geom_segment(aes(x=-60,y=0.4,xend= -60,yend=0.5),col='black',lwd=3)+
	geom_text(aes(x=-36,y=0.4),label='Novel altricials',size=16/.pt)+
	geom_text(aes(x=-46,y=0.4),label='Altricial',size=16/.pt)+
	geom_text(aes(x=-56,y=0.4),label='Precocial',size=16/.pt)+

	geom_segment(aes(x=-50,y=0.4,xend= -50,yend=0.5),col=cbp[2],lwd=3/2)+
	geom_segment(aes(x=-60,y=0.4,xend= -60,yend=0.5),col=cbp[3],lwd=3/2)+
	
	#geom_text_repel(data=df.sub, aes(y=ypos,x=0),angle=90, min.segment.length=1,parse=T, label=gsub('_.*','',df.sub$taxon),hjust=-1, nudge_x = .5, max.overlaps=3)+

	labs(col= expression(italic(sigma)), x='',y='')+
	geom_vline(xintercept=0.2,linewidth=2.5,colour='white')+
	theme(
	legend.position=c(0.10,0.20),
	legend.title=element_blank(),
         panel.background = element_rect(fill='white'),
         plot.background = element_rect(fill='white', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
		legend.background = element_rect(fill='transparent'),
		legend.key.width=unit( 0.8,'cm'),
		legend.key.height=unit( 1,'cm'),
		legend.text=element_text(size=16,color='black'),
		axis.line.x=element_line(color='black',linewidth=3/2),
		axis.line.y=element_line(color='black',linewidth=3/2),
		axis.ticks.x=element_line(color='black',linewidth=2),
		axis.ticks.y=element_line(color='black',linewidth=2),
		axis.title.x=element_text(size=30, color='black'),
		axis.title.y=element_text(size=30, color='black'),
		axis.text.x=element_text(size=16,color='black'),
		axis.text.y=element_text(size=16,color='black'),
)



# Repeat the above process for novel foot-uses instead of novel flight-styles 
dendr.bird.score <- cbind(dendr.trait,branch.col)

	dendr.bird.score <- cbind(dendr.bird.score,rep('NA',dim(dendr.trait)[1]))
	colnames(dendr.bird.score)[8]<-'nov'
	vector <- which(dendr$segments$yend==0)
	dendr.bird.score$nov[vector[which(lab.dat$label %in% Novel_feet==T)]]<-'Novel'

	dendr.novel<-dendr.bird.score
	dendr.novel$branch.col[which(dendr.bird.score$nov=='Novel')]<-'Novel Altricial'

df<-data.frame(Foot.use[,9:13])
rownames(df)<-Foot.use[,1]
df2<-cbind(rownames(df),df,dendr.bird.score$xend[which(dendr.bird.score$yend==0)][match(rownames(df), dendr$labels$label)] )
df2<-df2[complete.cases(df2),]
colnames(df2)[1]<-'taxon'
colnames(df2)[dim(df2)[2]]<-'ypos'
df2<-melt(df2, id=c('taxon','ypos'))
df2<-df2[-which(df2$value=='0'),]
df2$variable <- factor(df2$variable, levels= c(names(table(df2$variable))[order(table(df2$variable))]))
lab.order<-dendr$labels$label
df2$taxon <- factor(df2$taxon, levels= lab.order)
df2<-rbind(df2,df2)

df.sub <- df2[,1:2]
df.sub <-df.sub[!duplicated(df.sub), ]
df.sub <- df.sub[order(df.sub$ypos),]
df.sub <- df.sub[c(1,2,4,5,8,9,10,11,13,15,17,18,19,20),]



dev.tree.feet <- 
ggplot()+
	annotation_custom(grob=Upupa.g, xmin=0, xmax=20, ymin=-0.01214379-0.05, ymax=-0.01214379+0.05)+
	annotation_custom(grob=Topaza.g, xmin=0, xmax=20, ymin=-0.28959137-0.05, ymax=-0.28959137+0.05)+
	annotation_custom(grob=Pelecanus.g, xmin=2, xmax=22, ymin=0.32817671-0.05, ymax=0.32817671+0.05)+
	annotation_custom(grob=Pelagodroma.g, xmin=4, xmax=24, ymin=-0.39730291-0.05, ymax=-0.39730291+0.05)+
	annotation_custom(grob=Nyctibius.g, xmin=-6, xmax=16, ymin=0.26653139-0.02, ymax=0.26653139+0.04)+
	annotation_custom(grob=Leucocarbo.g, xmin=-5, xmax=20, ymin=0.10746715-0.01, ymax=0.10746715+0.1)+
	annotation_custom(grob=Jynx.g, xmin=0, xmax=20, ymin=-0.17065767-0.03, ymax=-0.17065767+0.03)+
	annotation_custom(grob=Fregata.g, xmin=0, xmax=20, ymin=0.50721715-0.08, ymax=0.50721715+0.05)+
	annotation_custom(grob=Columba.g, xmin=0, xmax=20, ymin=-0.12040668-0.02, ymax=-0.12040668+0.06)+
	xlim(-72.907238,19)+
	coord_flip()+
	geom_histogram(aes( y=ypos, fill=variable), data= df2, binwidth=0.02, alpha=.5)+
	scale_fill_manual(values= c( "#009E73", 
                      "#F0E442", "#0072B2", "#D55E00", "#CC79A7", 'darkgrey', 'lightgrey'),
labels=c('Raptorial','Scansorial','Digit manipulation','Forceful grip','Hanging clinging') )+
	geom_segment(data=dendr.bird.score[which(dendr.bird.score$branch.col=='Altricial'),],
  	aes(x = -y, y = x, xend = -yend, yend =xend ),col='black', lwd=3)+
	geom_segment(data=dendr.bird.score[which(dendr.bird.score$branch.col=='Altricial'),],
  	aes(x = -y, y = x, xend = -yend, yend =xend ),col=cbp[2], lwd=3/2)+

	geom_segment(data=dendr.bird.score[which(dendr.bird.score$branch.col=='Precocial'),],
  	aes(x = -y, y = x, xend = -yend, yend =xend ),col='black', lwd=3)+
	geom_segment(data=dendr.bird.score[which(dendr.bird.score$branch.col=='Precocial'),],
  	aes(x = -y, y = x, xend = -yend, yend =xend ),col=cbp[3], lwd=3/2)+

	geom_segment(data=dendr.bird.score[which(dendr.bird.score$nov=='Novel'),],
  	aes(x = -y, y = x, xend = -yend, yend =xend ),col='black', lwd=3/2)+

	#geom_text_repel(data=df.sub, aes(y=ypos,x=0),angle=90, min.segment.length=1,parse=T, label=gsub('_.*','',df.sub$taxon),hjust=-1, nudge_x = .5, max.overlaps=2)+

	#geom_text(data=df2, aes(y=ypos,x=7,label=taxon),angle=90)+
	#labs(col= expression(italic(sigma)), x=expression('Mya'), y=expression('(wing-leg)/(leg+wing)'))+
	labs(col= expression(italic(sigma)), x='',y='')+
	geom_vline(xintercept=0.2,linewidth=2.5,colour='white')+
	theme(
	legend.position=c(0.08,0.25),
	legend.title=element_blank(),
         panel.background = element_rect(fill='white'),
         plot.background = element_rect(fill='white', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
		legend.background = element_rect(fill='transparent'),
		legend.key.width=unit( 0.8,'cm'),
		legend.key.height=unit( 1,'cm'),
		legend.text=element_text(size=16,color='black'),
		axis.line.x=element_line(color='black',linewidth=3/2),
		axis.line.y=element_line(color='black',linewidth=3/2),
		axis.ticks.x=element_line(color='black',linewidth=2),
		axis.ticks.y=element_line(color='black',linewidth=2),
		axis.title.x=element_text(size=30, color='black'),
		axis.title.y=element_text(size=30, color='black'),
		axis.text.x=element_text(size=16,color='black'),
		axis.text.y=element_text(size=16,color='black'),
)

library(ggpubr) # 0.6.0

dev.new()
setwd('Z:/Andrew_backup/Precocial_paper')
prepped <- ggarrange(dev.tree.flight,dev.tree.feet,ncol=1,nrow=2, labels=c('a','b'), vjust=1, font.label = list(size = 30, color = "black", face = "bold", family = NULL))


labelled <- annotate_figure(prepped, left =  text_grob(expression("Mya"), 
              color = "black", face = "bold", size = 30, rot=90, vjust=.8))
labelled <- annotate_figure(labelled, bottom =text_grob(expression("(wing-leg)/(wing+leg)"), 
               color = "black", face = "bold", size = 30, vjust=-.75, hjust=0.3))

labelled
# 

ggsave(filename='FigIII_11_18_2024.pdf',height=45,width=45,unit='cm', dpi=300) 
# Save output 

# You can remove rows from df.sub manually later to make sure you onyl get the birds that you are interested in showing 

