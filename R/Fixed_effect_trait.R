library(metafor)
library(ggplot2)
library(ape)
library(rotl)
library(multcomp)
library(dplyr)
library(ggtree)

# To install the orchaRd package:
#install.packages("pacman")
#pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, emmeans)
#devtools::install_github("daniel1noble/orchaRd", force = TRUE)
library(orchaRd)
library(clubSandwich)
library(rmarkdown)

### Read in effect size data
effectdata <- read.csv(here("Data", "Survival project all pairwise.es.csv"))
effectdata <-  subset(effectdata, Paper.code != "HUM251") 


repdata_warm <- subset(effectdata, Trait.category == "Reproduction" & warm.cool == "Warm" )
repdata_cool <- subset(effectdata, Trait.category == "Reproduction" & warm.cool == "Cool" )

repdata <- rbind(repdata_warm, repdata_cool)



### Read in effect size data
longdata_warm <- subset(effectdata, Trait.category == "Longevity" & warm.cool == "Warm" )
longdata_cool <- subset(effectdata, Trait.category == "Longevity" & warm.cool == "Cool" )

alllong <- rbind(longdata_warm, longdata_cool)

rdata <- rbind(repdata, alllong)


########### change species names in survival data ####################################
classes <- read.csv("Data/Species_classifications.CSV") ## read in species classifications from map

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
### specify classifications from map
rdata$Class <- classes$class[match(rdata$Species.latin, classes$species_latin)]

### Create random factors into data frame 
rdata$obs <- factor(c(1:nrow(rdata)))                # Unique observation code
rdata$study_code <- factor(rdata$Paper.code)         # Model requires column names study_code (this is biggest level of nested code structure)
rdata$Species.phylo <- factor(rdata$Species.latin)   # Species names for phylo matrix
rdata$species <- factor(rdata$Species.latin)         # Another species column for ranom factor 

precision <- sqrt(1/rdata$v)                         # inverse standard error 
rdata[,"precision"] <- precision
str(rdata)

nlevels(rdata$species)    # Check number of species
nlevels(rdata$study_code) # Check number of studies


### setting covariance marix ###
rdata$shared_control <- factor(rdata$Effect.size.code)
VCV_shared <- impute_covariance_matrix(vi=rdata$v, cluster = rdata$shared_control, r=0.5)

#### running models ####

#### Subset main dataframe by experiments which have both longevity and reproduction
###  (these are the experiments in the unique_combo data)

uniquedata <- read.csv("Data/cleaned_unique_combo.csv")
expcodes <- as.vector(unique(uniquedata$Experiment.code))

test <- subset(rdata, Experiment.code %in% expcodes)

VCV_shared <- impute_covariance_matrix(vi=test$v, cluster = test$shared_control, r=0.5)



## Random factors: obs, study_code. Fixed Factors: Trait and treatment temp
meta_trait_test <- rma.mv(es, VCV_shared, mod=~Trait.category*treattemp,
                     random= list(~ 1|study_code, ~1|obs), data= test, method= "REML")
summary(meta_trait)







