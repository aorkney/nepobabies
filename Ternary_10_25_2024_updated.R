# The purpose of this script is to prouce ternary diagrams of fore and hindlimb proportions across 
# birds with distinct precocial and altricial developmental modes
# Much of this script is dedicated to formatting initial datasets. 
# I suggest we simply make a single spreadsheet of values for limb proportions and developmental mode
# when uploading this on GitHub. 
#
# Do you think it is acceptable for our GitHub deposition to include measurements derived from Bjarnason and Brinkworth?
# The code will be considerably longer if we have to derive all these features from Raw supplementary files from their publications. 
#
# There are penguins, ostriches, all sorts in the spaces

library(ggplot2) # version
# purpose
library(phytools)# version
# purpose
library(geomorph)# version
# purpose
library(nlme)# version
# purpose

# We want to plot both Bjarnason's and Brinkworth's unique bird species
# in a ternary diagram space for wings and legs, subdivided by altricial and precocial developmental strategy

setwd('C:/Users/Lab/Documents/AOrkney/bird_scripts_and_data_15_02_2023')
# Set work directory appropriately 

# Load requisite datasets
load('tree.22.10.2022.RData') # Load a phylogeny (we need this for technical reasons, such as matching congeners)
load('tree.names.22.10.2022.RData') # The names of the specimens from Bjarnason will be matched to congeners on the Prum tree for consistency with the rest of the study
load('coords.22.10.2022.RData') # Landmark constellations from Bjarnason et al., 2020
load('Csize.22.10.2022.RData') # We need these data to adjust the size of landmark constellations to their original volumes.
metadata <- read.csv("altriciality_scores_10_15_2023.2.csv") # Accompanying metadata for developmental strategy scores 



keeps <- dimnames(GPA.coords[[1]])[[3]]
# Original specimen names (not congener matches to the tree of Prum et al., 2015)

# Ensure names match the phylogeny.
# We are taking this step because the developmental modes in the metadata are associated with species labels found in the tree
# of Prum et al., 2015
for( i in 1:length(GPA.Csize) ){
  names(GPA.Csize[[i]]) <- tree_names 
	dimnames(GPA.coords[[i]])[[3]] <- tree_names
}


dfCsize.birds <- do.call(cbind,GPA.Csize)
# Format the data for centroid sizes into a matrix

developmental_mode <- metadata$my[match(tree_names,metadata$tree)]
developmental_mode[grep('emi', developmental_mode)]<-'Semi'
developmental_mode[grep('tricial', developmental_mode)]<-'Altricial'
developmental_mode[grep('cocial', developmental_mode)]<-'Precocial'
names(developmental_mode)<-tree_names
# Define developmental strategy for each species

tips <- tree_names[which(developmental_mode == 'Precocial' | developmental_mode == 'Altricial')]
# Select species names from the tree of Prum et al., 2015, for which developmental mode is explicitly precocial or altricial

# Define a fuction to find euclidean distances between points
euclid<-function(data){
  return(sqrt( (data[1,1]-data[2,1])^2 + (data[1,2]-data[2,2])^2 + (data[1,3]-data[2,3])^2 ))
}

# Define a function that can find the euclidean distances between defined points in a landmark constellaiton
# representing a bone; we will need this in order to compute bone lengths from the Bjarnason dataset. 
find.euclid.distance <- function(bone,coords,Csize,landmarks){
  distances<-list()
  for(i in 1:dim(coords[[bone]])[3]){
    inflated <- coords[[bone]][,,i]*Csize[[bone]][i]
    distances[[i]]<-euclid(inflated[landmarks,])
  }
  distances<-unlist(distances)
  names(distances)<-names(Csize[[bone]])
  return(distances)
}

humerus.length.Bj <- find.euclid.distance('humerus',GPA.coords,GPA.Csize,landmarks=c(4,70)) 
radius.length.Bj <- find.euclid.distance('radius',GPA.coords,GPA.Csize,landmarks=c(9,3)) 
ulna.length.Bj <- find.euclid.distance('ulna',GPA.coords,GPA.Csize,landmarks=c(8,3)) 
carpometacarpus.length.Bj <- find.euclid.distance('carpometacarpus',GPA.coords,GPA.Csize,landmarks=c(16,70)) 
femur.length.Bj <- find.euclid.distance('femur',GPA.coords,GPA.Csize,landmarks=c(42,23))
tibiotarsus.length.Bj <- find.euclid.distance('tibiotarsus',GPA.coords,GPA.Csize,landmarks=c(15,49))
tarsometatarsus.length.Bj <- find.euclid.distance('tarsometatarsus',GPA.coords,GPA.Csize,landmarks=c(1,50))
# Generate bone length measurements from landmark constellations; the individual landmark numbers
# were chosen after visual inspection to affirm they represent valid length measurements. 

developmental_mode.Bj<-developmental_mode
# Store developmental strategy information in an object associated with the Bjarnason dataset. 

names(humerus.length.Bj) <- keeps 
names(radius.length.Bj) <- keeps 
names(ulna.length.Bj) <- keeps 
names(carpometacarpus.length.Bj) <- keeps 
names(femur.length.Bj) <- keeps 
names(tibiotarsus.length.Bj) <- keeps 
names(tarsometatarsus.length.Bj) <- keeps 
names(developmental_mode.Bj) <- keeps
# Restore original specimen species names. 

# Thus far we have used only Bjarnason 2021's dataset, 
# but we can also combine this with the Brinkworth et al., 2023 aggregated suite of data. 

# Load Brinkworth's aggregated dataset. 
setwd('C:/Users/Lab/Documents/AOrkney/Precocial_study')
brinkworth.data <- read.csv('altriciality_scores_brinkworth_06_03_2024.csv')
# Load the brinkworth data and its accompanying developmental strategy scores. 

setwd('C:/Users/Lab/Documents/AOrkney/Brinkworth/23941488/Supplementary_Software/Supplementary_Software')
jetz.tree <- read.tree('Prum_merge_hackett_stage2_1K_mcc.tree')
# Load a corresponding tree for the Brinkworth dataset; we will use Jetz et al., 2010

pruned.tree <- keep.tip(jetz.tree, brinkworth.data$JetzTreeName)
# Prune the tree to the taxa of interest; this requires some matching of specimens in Brinkworth to congeners in the Jetz tree

GPA.Csize <- list()
GPA.Csize[[1]] <- brinkworth.data[,5]
GPA.Csize[[2]] <- brinkworth.data[,7]
GPA.Csize[[3]] <- brinkworth.data[,8]
GPA.Csize[[4]] <- brinkworth.data[,9]
GPA.Csize[[5]] <- brinkworth.data[,10]
GPA.Csize[[6]] <- brinkworth.data[,11]

for(i in 1:length(GPA.Csize)){
	names(GPA.Csize[[i]]) <- brinkworth.data$Jetz
}
# Compile a matrix of limb length measurements from the Brinkworth dataset

names(GPA.Csize) <- c('humerus','radius','carpometacarpus','femur','tibiotarsus','tarsometatarsus')
masses <- brinkworth.data$Body
names(masses)<-brinkworth.data$Jetz
# Render names consistent

developmental_mode <- brinkworth.data$altr[match(names(masses),brinkworth.data$Jetz)]
developmental_mode[grep('emi', developmental_mode)]<-'Semi'
developmental_mode[grep('tricial', developmental_mode)]<-'Altricial'
developmental_mode[grep('cocial', developmental_mode)]<-'Precocial'
names(developmental_mode)<-names(masses)
# Extract the developmental strategy categorisations and standardise them.  

developmental_mode.Br <- developmental_mode
# Specifically store the vector of developmental strategies in a distinct object associated with the Brinkworth dataset

humerus.length.Br <- GPA.Csize$humerus
radius.length.Br <- GPA.Csize$radius
ulna.length.Br <- GPA.Csize$ulna
carpometacarpus.length.Br <- GPA.Csize$carpometacarpus
femur.length.Br <- GPA.Csize$femur
tibiotarsus.length.Br <- GPA.Csize$tibiotarsus
tarsometatarsus.length.Br <- GPA.Csize$tarsometatarsus
# Bone lengths derived from the Brinkworth dataset. 
# No name matches are necessary here because aall of the original specimen names from the Brinkworth dataset were represented in the 
# Jetz tree. 

overlap <- names(humerus.length.Br) [which(names(humerus.length.Br) %in% names(humerus.length.Bj) ==T)]
# Compute overlap between the species in the Brinkworth and Bjarnason dataset

plot(tarsometatarsus.length.Br[overlap],
tarsometatarsus.length.Bj[overlap])
abline(0,1)
# Assess whether datasets have similar values for the overlapping species


humerus.length <- c(humerus.length.Bj,humerus.length.Br[which(names(humerus.length.Br) %in% overlap==F)])
radius.length <- c(radius.length.Bj,radius.length.Br[which(names(humerus.length.Br) %in% overlap==F)])
carpometacarpus.length <- c(carpometacarpus.length.Bj,carpometacarpus.length.Br[which(names(humerus.length.Br) %in% overlap==F)])
femur.length <- c(femur.length.Bj,femur.length.Br[which(names(humerus.length.Br) %in% overlap==F)])
tibiotarsus.length <- c(tibiotarsus.length.Bj,tibiotarsus.length.Br[which(names(humerus.length.Br) %in% overlap==F)])
tarsometatarsus.length <- c(tarsometatarsus.length.Bj,tarsometatarsus.length.Br[which(names(humerus.length.Br) %in% overlap==F)])
# Pool the Bjarnason and Brinkworth data and eliminate overlapping taxa from the Brinkworth arm of the dataset. 

df<-data.frame(cbind(humerus.length,radius.length,carpometacarpus.length,femur.length,tibiotarsus.length,tarsometatarsus.length))
colnames(df)<-c('humerus','radius','carpometacarpus','femur','tibiotarsus','tarsometatarsus')
df[,1:3]<-df[,1:3]/rowSums(df[,1:3])
df[,4:6]<-df[,4:6]/rowSums(df[,4:6])
df<-cbind(df,c(developmental_mode.Bj,developmental_mode.Br[which(names(humerus.length.Br) %in% overlap==F)]) )
colnames(df)[7]<-'mode'
# Bind the data into a dataframe and add developmental strategy categorisations. 


# The data are prepared for plotting. The following code will produce visualisations. 

# Function to convert degrees to radians 
deg2rad <- function(deg) {(deg * pi) / (180)}
tern2cart <- function(m1,m2){
	x <- -(m1 + (m2)*cos(deg2rad(60)))
	y <- m2*sin(deg2rad(60))
	return(cbind(x,y))
}

segdat <- data.frame(tern2cart(c( 100,0,  0,0, 0,100 ),c( 0,100, 100,0, 0,0 )))
outline <- data.frame(cbind(segdat$x[seq(1,dim(segdat)[1],by=2)],segdat$x[seq(2,dim(segdat)[1],by=2)],segdat$y[seq(1,dim(segdat)[1],by=2)],segdat$y[seq(2,dim(segdat)[1],by=2)]))
segdat <- data.frame(tern2cart(c(55,35, 35,15, 15,15, 15,35, 35,55, 55,55  ),c(35,55, 55,55, 55,35, 35,15, 15,15, 15,35  )))
segments <- data.frame(cbind(segdat$x[seq(1,dim(segdat)[1],by=2)],segdat$x[seq(2,dim(segdat)[1],by=2)],segdat$y[seq(1,dim(segdat)[1],by=2)],segdat$y[seq(2,dim(segdat)[1],by=2)]))
# Make a framework of lines to outline the ternary plot 

#The following code makes a dense network of lines spaced at 10% increments throughout the ternary plot: 
x <- rep(seq(90,10,length.out=9),each=2)
x0 <- x
x0[seq(1,length(x),2)]<-0
xS <-rep(seq(10,90,length.out=9),each=2)
xS[seq(1,length(x),2)]<-0
x <- c(x, x0, xS)
y <- rep(seq(10,90,length.out=9),each=2)
y[seq(1,length(y),2)]<-0
yS<-rep(seq(10,90,length.out=9),each=2)
yS[seq(2,length(yS),2)]<-0
y <- c(y, rep(seq(10,90,length.out=9),each=2),yS)


# I suspect Argusianus_argus bird 555, from Brinkworth's data, has been measured incorrectly; its carpometacarpus is unusually large. 
which(rownames(df)=='Argusianus_argus')
df<-df[-556,]
# Remove this bird. 

# The following code will produce lines and tick marks to partition the ternary diagram. 

segdat <- data.frame(tern2cart(x,y))
colnames(segdat)<-c('x','y')
lines <- data.frame(cbind(segdat$x[seq(1,dim(segdat)[1],by=2)],segdat$x[seq(2,dim(segdat)[1],by=2)],segdat$y[seq(1,dim(segdat)[1],by=2)],segdat$y[seq(2,dim(segdat)[1],by=2)]))
# Produce lines that divide the ternary diagram into deciles.

# Ticks along base axis:
tickdat <- data.frame(tern2cart(rep(seq(90,10,length.out=9),each=2),rep(0,18) ))
colnames(tickdat) <- c('x','y')
ticks <- data.frame(cbind(tickdat$x[seq(1,dim(tickdat)[1],by=2)],tickdat$x[seq(2,dim(tickdat)[1],by=2)],
tickdat$y[seq(1,dim(tickdat)[1],by=2)],tickdat$y[seq(2,dim(tickdat)[1],by=2)] ) )
ticks$X4 <- ticks$X4 -2
ticks$X2 <- ticks$X2 -2*cos(deg2rad(60))

# Ticks along right hand slope of ternary plot:
ticks2 <- ticks
ticks2$X1 <- seq(-10,-90,length.out=9)*cos(deg2rad(60))
ticks2$X2 <- ticks2$X1 +2
ticks2$X3 <- seq(10,90,length.out=9)*sin(deg2rad(60))
ticks2$X4 <- seq(10,90,length.out=9)*sin(deg2rad(60))

# Ticks along left hand slope of ternary plot:
ticks3 <- ticks
ticks3$X1 <- seq(-90,-10,length.out=9)*cos(deg2rad(60))-50
ticks3$X2 <- ticks3$X1-cos(deg2rad(60))*2
ticks3$X3 <- seq(10,90,length.out=9)*sin(deg2rad(60))
ticks3$X4 <- seq(10,90,length.out=9)*sin(deg2rad(60))+sin(deg2rad(60))*2

# Arrows that extend around the border of the plot to indicate axial directions:
Arrows <- data.frame( c(0,0,0), c(0,0,0), c(0,0,0), c(0,0,0) )
colnames(Arrows) <- c('X1','X2','X3','X4')
Arrows$X1[1] <- 4
Arrows$X3[1] <- 0
Arrows$X2[1] <- 4 - 20*cos(deg2rad(60))
Arrows$X4[1] <- 20 *sin(deg2rad(60))
Arrows$X1[2] <- -50 - 4*cos(deg2rad(60))
Arrows$X3[2] <- 90
Arrows$X2[2] <- -70 *sin(deg2rad(60)) - 4*cos(deg2rad(60))
Arrows$X4[2] <- 90 - (20 *sin(deg2rad(60)))
Arrows$X1[3] <- -100 -4*cos(deg2rad(60))
Arrows$X3[3] <- -2*tan(deg2rad(60))
Arrows$X2[3] <- -80
Arrows$X4[3] <- -2*tan(deg2rad(60))

# Produce a ternary diagram, with a shaded hexagonal region containing all the flying birds'
# forelimb proportions. (two examples of flightless birds lie beyond this space; the ostrich and the cassowary)
# A network of fine hexagons illustrate the occupied region of this space. 

wing.base <- 
ggplot( ) +
geom_segment(data=ticks,aes(x=X1,xend=X2,y=X3,yend=X4),linewidth=1,col='black',alpha=1)+
geom_segment(data=ticks2,aes(x=X1,xend=X2,y=X3,yend=X4),linewidth=1,col='black',alpha=1)+
geom_segment(data=ticks3,aes(x=X1,xend=X2,y=X3,yend=X4),linewidth=1,col='black',alpha=1)+
geom_segment(data=Arrows,aes(x=X1,xend=X2,y=X3,yend=X4),linewidth=1,col='black',alpha=1,  arrow = arrow(length = unit(0.5, "cm")), lineend='round',linejoin='round')+
geom_segment(data=lines,aes(x=X1,xend=X2,y=X3,yend=X4),linewidth=1/2,col='lightgrey',alpha=0.5)+
geom_segment(data=outline,aes(x=X1,xend=X2,y=X3,yend=X4),linewidth=3/4,col='black')+
geom_hex(data=data.frame(tern2cart(df$humerus*100,df$radius*100),df$mode),aes(x=x,y=y),binwidth=1,fill=NA,col='black',linewidth=1/4)+
#geom_text(data=data.frame(tern2cart(df$humerus*100,df$radius*100),df$mode),aes(x=x,y=y), label=rownames(df))+
geom_polygon(data=segments,aes(x=X1,y=X3),linewidth=3/4,col='black',alpha=0.2)+
geom_text( data=data.frame(cbind(c(-102-2,-44+2,-10),c(3+2,82.60254+2,-2-2))) , aes(x=X1,y=X2, label=c('Humerus','Radius','Carpometacarpus')),angle=c(60,-60,0),size=8/.pt)+ 
theme(  panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	axis.ticks.x=element_blank(),
	legend.position='none' )+
 coord_fixed()


# Produce a plot for altricial birds within the shaded hexagon, with a shaded density distribution
segdat <- data.frame(tern2cart(c(55,35, 35,15, 15,15, 15,35, 35,55, 55,55  ),c(35,55, 55,55, 55,35, 35,15, 15,15, 15,35  )))
segments <- data.frame(cbind(segdat$x[seq(1,dim(segdat)[1],by=2)],segdat$x[seq(2,dim(segdat)[1],by=2)],segdat$y[seq(1,dim(segdat)[1],by=2)],segdat$y[seq(2,dim(segdat)[1],by=2)]))
dfalt <-data.frame(tern2cart(df$humerus*100,df$radius*100))[which(df$mode=='Altricial'),]

# You need to crop out the two ratite birds?
hexies <- data.frame(tern2cart(df$humerus*100,df$radius*100),df$mode)
hexies <- hexies[-which(df$humerus > .55),]

wing.alt <- 
ggplot( ) +
  stat_density_2d(data=dfalt, aes(x=x,y=y,fill = ..ndensity.. ,  alpha = ..ndensity..^.5), geom = "raster", contour = F)+
scale_fill_continuous(limits = c(0, 1), type = "viridis",option='plasma')+
scale_alpha_continuous(range = c(0.0, 1), limits = c(0, 0.5), guide = guide_none()) +
geom_segment(data=segments,aes(x=X1,xend=X2,y=X3,yend=X4),linewidth=3/2,col='#E69F00',lineend="round")+
#lims(x=c(-80,-30))+
geom_hex(data=data.frame(hexies),aes(x=x,y=y),binwidth=1,fill=NA,col='black',linewidth=1/4)+
theme(  panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	axis.ticks.x=element_blank(),
	legend.position='none' )+
 coord_fixed()


# Produce a plot for semi-altricial or semi-precocial birds within the shaded hexagon, with a shaded density distribution
# We won't actually present this in the final compiled figure but some readers may be interested
dfsemi <-data.frame(tern2cart(df$humerus*100,df$radius*100))[grep('emi',df$mode),]
wing.semi <- 
ggplot( ) +
  stat_density_2d(data=dfsemi, aes(x=x,y=y,fill = ..ndensity.. ,  alpha = ..ndensity..^.5), geom = "raster", contour = F)+
scale_fill_continuous(limits = c(0, 1), type = "viridis",option='plasma')+
scale_alpha_continuous(range = c(0.0, 1), limits = c(0, 0.5), guide = guide_none()) +
geom_segment(data=segments,aes(x=X1,xend=X2,y=X3,yend=X4),linewidth=3/2,col='black',lineend="round")+
#lims(x=c(-80,-30))+
geom_hex(data=data.frame(hexies),aes(x=x,y=y),binwidth=1,fill=NA,col='black',linewidth=1/4)+
theme(  panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	axis.ticks.x=element_blank(),
	legend.position='none' )+
 coord_fixed()

# Produce a plot for precocial birds within the shaded hexagon, with a shaded density distribution
dfpre <-data.frame(tern2cart(df$humerus*100,df$radius*100))[grep('recocial',df$mode),]
dfpre <- dfpre[-which(df[grep('recocial',df$mode),]$humerus >0.55),]
wing.pre <- 
ggplot( ) +
  stat_density_2d(data=dfpre, aes(x=x,y=y,fill = ..ndensity.. ,  alpha = ..ndensity..^.5), geom = "raster", contour = F)+
scale_fill_continuous(limits = c(0, 1), type = "viridis",option='plasma')+
scale_alpha_continuous(range = c(0.0, 1), limits = c(0, 0.5), guide = guide_none()) +
geom_segment(data=segments,aes(x=X1,xend=X2,y=X3,yend=X4),linewidth=3/2,col='#56B4E9',lineend='round')+
#lims(x=c(-80,-30))+
geom_hex(data=data.frame(hexies),aes(x=x,y=y),binwidth=1,fill=NA,col='black',linewidth=1/4)+
theme(  panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	axis.ticks.x=element_blank(),
	legend.position='none' )+
 coord_fixed()


library(cowplot)
# We need this package in order to arrange multiple subplots

# Produce a blank plot we can use as a backdrop
blank<-ggplot()+ geom_blank()+
  theme(
	axis.text.x=element_text(size=7,colour='black'),
  	axis.text.y=element_text(size=7,colour='black'),
	axis.title.x.bottom=element_text(size=7,colour="black"),
	axis.title.y.left=element_text(size=7,colour="black"),
      axis.ticks=element_line(size=1),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())


# Begin constructing a baseplot for the leg ternary diagram
segdat <- data.frame(tern2cart(c( 100,0,  0,0, 0,100 ),c( 0,100, 100,0, 0,0 )))
outline <- data.frame(cbind(segdat$x[seq(1,dim(segdat)[1],by=2)],segdat$x[seq(2,dim(segdat)[1],by=2)],segdat$y[seq(1,dim(segdat)[1],by=2)],segdat$y[seq(2,dim(segdat)[1],by=2)]))
segdat <- data.frame(tern2cart(c( 45,45, 45,30, 30,5, 5,5, 5,20, 20,45  ),c( 30,45, 45,60, 60,60, 60,45, 45,30, 30,30 )))
segments <- data.frame(cbind(segdat$x[seq(1,dim(segdat)[1],by=2)],segdat$x[seq(2,dim(segdat)[1],by=2)],segdat$y[seq(1,dim(segdat)[1],by=2)],segdat$y[seq(2,dim(segdat)[1],by=2)]))

x <- rep(seq(90,10,length.out=9),each=2)
x0 <- x
x0[seq(1,length(x),2)]<-0
xS <-rep(seq(10,90,length.out=9),each=2)
xS[seq(1,length(x),2)]<-0
x <- c(x, x0, xS)
y <- rep(seq(10,90,length.out=9),each=2)
y[seq(1,length(y),2)]<-0
yS<-rep(seq(10,90,length.out=9),each=2)
yS[seq(2,length(yS),2)]<-0
y <- c(y, rep(seq(10,90,length.out=9),each=2),yS)

segdat <- data.frame(tern2cart(x,y))
colnames(segdat)<-c('x','y')
lines <- data.frame(cbind(segdat$x[seq(1,dim(segdat)[1],by=2)],segdat$x[seq(2,dim(segdat)[1],by=2)],segdat$y[seq(1,dim(segdat)[1],by=2)],segdat$y[seq(2,dim(segdat)[1],by=2)]))

hexies <- data.frame(tern2cart(df$femur*100,df$tibiotarsus*100),df$mode)
hexies <- hexies[-which(df$femur > .40),]
# Remove Loriculus galgulus because I think it was measured badly 

# Construct the leg proportion ternary diagram
leg.base <-
ggplot( ) +
geom_segment(data=ticks,aes(x=X1,xend=X2,y=-X3,yend=-X4),linewidth=1,col='black',alpha=1)+
geom_segment(data=ticks2,aes(x=X1,xend=X2,y=-X3,yend=-X4),linewidth=1,col='black',alpha=1)+
geom_segment(data=ticks3,aes(x=X1,xend=X2,y=-X3,yend=-X4),linewidth=1,col='black',alpha=1)+
geom_segment(data=Arrows,aes(x=X1,xend=X2,y=-X3,yend=-X4),linewidth=1,col='black',alpha=1,  arrow = arrow(length = unit(0.5, "cm")), lineend='round',linejoin='round')+
geom_segment(data=lines,aes(x=X1,xend=X2,y=-X3,yend=-X4),linewidth=1/2,col='lightgrey',alpha=0.5)+
geom_segment(data=outline,aes(x=X1,xend=X2,y=-X3,yend=-X4),linewidth=3/4,col='black')+
geom_hex(data=hexies,aes(x=x,y=-y),binwidth=1,fill=NA,col='black',linewidth=1/4)+
geom_polygon(data=segments,aes(x=X1,y=-X3),linewidth=3/4,col='black',alpha=0.2)+
geom_text( data=data.frame(cbind(c(-100-2,-45+2,-8),c(3,82.60254,-2-2))) , aes(x=X1,y=-X2, label=c('Femur','Tibiotarsus','Tarsometatarsus')),angle=c(-60,60,0),size=8/.pt)+ 
#geom_text(data=data.frame(tern2cart(df$femur*100,df$tibiotarsus*100),df$mode)[grep('Apteno',rownames(df)),],aes(x=x,y=-y), label=rownames(df)[grep('Apteno',rownames(df))])+
theme(  panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	axis.ticks.x=element_blank(),
	legend.position='none')+
 coord_fixed()


segdat <- data.frame(tern2cart(c( 45,45, 45,30, 30,5, 5,5, 5,20, 20,45  ),c( 30,45, 45,60, 60,60, 60,45, 45,30, 30,30 )))
segments <- data.frame(cbind(segdat$x[seq(1,dim(segdat)[1],by=2)],segdat$x[seq(2,dim(segdat)[1],by=2)],segdat$y[seq(1,dim(segdat)[1],by=2)],segdat$y[seq(2,dim(segdat)[1],by=2)]))

# Construct a ternary diagram for the distribution of precocial species
dfpre <-data.frame(tern2cart(df$femur*100,df$tibiotarsus*100))[grep('ocial',df$mode),]
leg.pre<-
ggplot( ) +
  stat_density_2d(data=dfpre, aes(x=x,y=-y,fill = ..ndensity.. ,  alpha = ..ndensity..^.5), geom = "raster", contour = F)+
scale_fill_continuous(limits = c(0, 1), type = "viridis",option='plasma')+
scale_alpha_continuous(range = c(0.0, 1), limits = c(0, 0.5), guide = guide_none()) +
geom_segment(data=segments,aes(x=X1,xend=X2,y=-X3,yend=-X4),linewidth=3/2,col='#56B4E9',lineend='round')+
#lims(x=c(-70,-25),y=c(25,55))+
geom_hex(data=hexies,aes(x=x,y=-y),binwidth=1,fill=NA,col='black',linewidth=1/4)+
theme(  panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	axis.ticks.x=element_blank(),
	legend.position='none' )+
 coord_fixed()


# Construct a ternary diagram for the distribution of semialtricial or semiprecocial species
# We won't actually use this for the final plot but some readers might be interested. 
dfsemi <-data.frame(tern2cart(df$femur*100,df$tibiotarsus*100))[grep('emi',df$mode),]
leg.semi<-
ggplot( ) +
  stat_density_2d(data=dfsemi, aes(x=x,y=-y,fill = ..ndensity.. ,  alpha = ..ndensity..^.5), geom = "raster", contour = F)+
scale_fill_continuous(limits = c(0, 1), type = "viridis",option='plasma')+
scale_alpha_continuous(range = c(0.0, 1), limits = c(0, 0.5), guide = guide_none()) +
geom_segment(data=segments,aes(x=X1,xend=X2,y=-X3,yend=-X4),linewidth=3/2,col='#56B4E9',lineend='round')+
#lims(x=c(-70,-25),y=c(25,55))+
geom_hex(data=hexies,aes(x=x,y=-y),binwidth=1,fill=NA,col='black',linewidth=1/4)+
theme(  panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	axis.ticks.x=element_blank(),
	legend.position='none' )+
 coord_fixed()

# Construct a ternary diagram for the distribution of altricial species
dfalt <-data.frame(tern2cart(df$femur*100,df$tibiotarsus*100))[grep('ricial',df$mode),]
leg.alt<-
ggplot( ) +
  stat_density_2d(data=dfalt, aes(x=x,y=-y,fill = ..ndensity.. ,  alpha = ..ndensity..^.5), geom = "raster", contour = F)+
scale_fill_continuous(limits = c(0, 1), type = "viridis",option='plasma')+
scale_alpha_continuous(range = c(0.0, 1), limits = c(0, 0.5), guide = guide_none()) +
geom_segment(data=segments,aes(x=X1,xend=X2,y=-X3,yend=-X4),linewidth=3/2,col='#E69F00',lineend='round')+
lims(x=range(c(segments$X1,segments$X2)),y=range(c(-segments$X3,-segments$X4)))+
geom_hex(data=hexies,aes(x=x,y=-y),binwidth=1,fill=NA,col='black',linewidth=1/4)+
theme(  panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	axis.ticks.x=element_blank(),
	legend.position='none')+
	coord_fixed()


dev.new(width=12,height=10,unit='cm')
n <- 40/25.98076
# Selecting a suitable n will fix the aspect ratio 
ggdraw()+
draw_plot(blank)+
draw_plot(leg.base,x=0,y=0,height=0.5,width=0.5)+
draw_plot(wing.base,x=0,y=0.5,height=0.5,width=0.5)+
draw_plot(wing.pre, x=0.35,y=0.65,width=0.35,height=0.35)+
draw_plot(wing.alt, x=0.6,y=0.5,width=0.35,height=0.35)+
draw_plot(leg.pre, x=0.01,y=0.0,height=0.35*0.8)+
draw_plot(leg.alt, x= (0.35*0.8) - 0.005,y=0.2,height=0.35*0.8)+
# You also want labels for subplots
geom_text(aes(x=c(0.05,0.42,0.7,0.05,0.42,0.65),y=c(0.9,0.95,0.85,0.1,0.285,0.475)),label=c('a','b','c','d','e','f'),size=20/.pt,fontface='bold')
setwd('C:/Users/Lab/Documents/AOrkney/Precocial_study')
ggsave(filename='Ternaries_test_10_29_2024b.pdf',width=18,height=18,unit='cm') 

# Plot has been arranged and saved 
