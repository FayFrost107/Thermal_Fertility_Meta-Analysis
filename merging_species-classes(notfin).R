###-----------------------------------------------------###
############################################ Preamble ######################################################
rm(list=ls()) # Clear R environment

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

#--------------#
# 1. Setup     #
#--------------#

### Read in effect size data
effectdata <- read.csv("Data/Survival project all pairwise.es.csv")
longdata_warm <- subset(effectdata, Trait.category == "Longevity" & warm.cool == "Warm" )
longdata_cool <- subset(effectdata, Trait.category == "Longevity" & warm.cool == "Cool" )

alllong <- rbind(longdata_warm, longdata_cool)

### select data for analysis
rdata <- alllong

rdata <- subset(rdata, Paper.code != "HUM251")

rdata <- rdata %>% mutate(c_treattemp = treattemp - 25)


########### change species names in survival data ####################################
classes <- read.csv("Data/Systematic map species list.CSV") ## read in species classifications from map

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
classes$species_latin <- gsub("_", " ", classes$species_latin)
rdata$Class <- classes$class[match(rdata$Species.latin, classes$species_latin)]
rdata$Fertilisation.mode <- classes$fertilisation_mode[match(rdata$Species.latin, classes$species_latin)]
















