
rm(list=ls())

library(rotl)
library(tidyr)
library(dplyr)
library(ape)
library(ggplot2)
library(ggtree)


############################# import full map tree ################################################
tree1 <- read.nexus("Phylogeny/complete_tree_draft.nex")
tree_grafen = compute.brlen(tree1, method="Grafen", power=1)


### Read in effect size data
effectdata <- read.csv("Data/Survival project all pairwise.es.csv")
repdata_warm <- subset(effectdata, Trait.category == "Reproduction" & warm.cool == "Warm" )
repdata_cool <- subset(effectdata, Trait.category == "Reproduction" & warm.cool == "Cool" )

longdata_warm <- subset(effectdata, Trait.category == "Longevity" & warm.cool == "Warm" )
longdata_cool <- subset(effectdata, Trait.category == "Longevity" & warm.cool == "Cool" )

survdata_warm <- subset(effectdata, Trait.category == "Survival" & warm.cool == "Warm" )
survdata_cool <- subset(effectdata, Trait.category == "Survival" & warm.cool == "Cool" )

allrep <- rbind(repdata_warm, repdata_cool)
alllong <- rbind(longdata_warm, longdata_cool)
allsurv <- rbind(survdata_warm, survdata_cool)

### select data for tree
rdata <- allsurv
rdata <- subset(rdata, Paper.code != "HUM251")
rdata <- subset(rdata, Paper.code != "OSM205")

#rdata <- subset(effectdata, warm.cool != "Reference")   <<<<< select for all data
########### change species names in survival data

rdata$Species.latin[which(rdata$Species.latin == "Marasmia exigua")]                <- "Cnaphalocrocis exigua"
rdata$Species.latin[which(rdata$Species.latin == "Matsumuratettix hieroglyphicus")] <- "Matsumuratettix hiroglyphicus"
rdata$Species.latin[which(rdata$Species.latin == "Mythimna roseilinea")]            <- "Mythimna albipuncta"
rdata$Species.latin[which(rdata$Species.latin == "Apis craccivora")]                <- "Aphis craccivora"
rdata$Species.latin[which(rdata$Species.latin == "Cryptoleamus montrouzieri")]      <- "Cryptolaemus montrouzieri"
rdata$Species.latin[which(rdata$Species.latin == "Asplanchna brightwelli")]         <- "Asplanchna brightwellii"
rdata$Species.latin[which(rdata$Species.latin == "Brennandania lambi")]             <- "Pygmephorus lambi"
rdata$Species.latin[which(rdata$Species.latin == "Amblyseius alstoniae")]           <- "Euseius alstoniae"
rdata$Species.latin[which(rdata$Species.latin == "Siphoninus phyllyreae")]          <- "Siphoninus phillyreae"
rdata$Species.latin[which(rdata$Species.latin == "Proprioseiopsis asetus")]         <- "Amblyseius asetus"
rdata$Species.latin[which(rdata$Species.latin == "Parabemisia myrica")]             <- "Parabemisia myricae"
rdata$Species.latin[which(rdata$Species.latin == "Cirrospilus sp. near lyncus")]    <- "Cirrospilus lyncus"
rdata$Species.latin[which(rdata$Species.latin == "Anagyrus sp. nov. nr. sinope" )]  <- "Anagyrus sinope"
rdata$Species.latin[which(rdata$Species.latin == "Monochamus leuconotus")]          <- "Anthores leuconotus"
rdata$Species.latin[which(rdata$Species.latin == "Ropalosiphum maidis")]            <- "Rhopalosiphum maidis"

rdata$Species.latin[which(rdata$Species.latin == "Artemia fransiscana")]            <- "Artemia franciscana"
rdata$Species.latin[which(rdata$Species.latin == "Blathyplectes curculionis")]      <- "Bathyplectes curculionis"
rdata$Species.latin[which(rdata$Species.latin == "Menochilus sexmaculatus")]        <- "Cheilomenes sexmaculata"
rdata$Species.latin[which(rdata$Species.latin == "unknown (Tominic)")]              <- "Trichogramma" 




### Problem species from systematic map to classes
rdata$Species.latin[which(rdata$Species.latin == "Cosmocomoidea ashmeadi")]              <- 	"Gonatocerus ashmeadi"
rdata$Species.latin[which(rdata$Species.latin == "Cosmocomoidea triguttata")]   <- "Gonatocerus triguttatus"
rdata$Species.latin[which(rdata$Species.latin == "Mythimna albipuncta")]            <-  "Mythimna roseilinea"
rdata$Species.latin[which(rdata$Species.latin == "Daphnia australis")]            <-  "Daphniopsis australis"

## import tree from map
notin <- setdiff(tree1$tip.label, rdata$Species.latin)
indata <- setdiff(tree1$tip.label, notin)

## number of indata should equal # of unique rdata species. 
setdiff(unique(rdata$Species.latin), indata)


## prune tree
prune_tree <- drop.tip(tree_grafen, notin)



write.nexus(prune_tree, file="complete_tree.nex")


##### create data frame and plot tree ###################
classes <- read.csv("Species_classifications.CSV") ## read in species classifications from map

plot_data <- c()
plot_data$species_latin <- indata
plot_data$class <- rdata$Class[match(indata, rdata$Species.latin)]
plot_data <- as.data.frame(unique(plot_data))
colnames(plot_data) <- c("species_latin", "class")




ggtree(prune_tree, layout = "circular", lwd = 0.1) %<+% plot_data +
  geom_tiplab(size=2, offset=0.01) +
  aes(col = class) 
