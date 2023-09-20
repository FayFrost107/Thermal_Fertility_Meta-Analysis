###-----------------------------------------------------###
### Multilevel meta-analysis using Metafor              ###
### Author: Fay Frost [fay.frost@liverpool.ac.uk]             ###
### Code adapted from Liam Dougherty. 
### University of Liverpool                             ###
### Date: August 2023                                   ###
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

##########################################################################################################


#--------------#
# 1. Setup     #
#--------------#

### Read in effect size data
effectdata <- read.csv("./data/Survival project all pairwise.es.csv")
longdata_warm <- subset(effectdata, Trait.category == "Longevity" & warm.cool == "Warm" )
longdata_cool <- subset(effectdata, Trait.category == "Longevity" & warm.cool == "Cool" )

alllong <- rbind(longdata_warm, longdata_cool)

### select data for analysis
rdata <- alllong

rdata <- subset(rdata, Paper.code != "HUM251")


########### change species names in survival data ####################################
classes <- read.csv("./data/Species_classifications.CSV") ## read in species classifications from map

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

#### Import Tree #############

## import tree from map
tree1 <- read.nexus("./phylogeny/all_long_excHUM251_tree.nex")
tree_grafen = compute.brlen(tree1, method="Grafen", power=1)
phylo_matrix <- vcv(tree_grafen, cor=TRUE, model="Brownian") # Make phylogenetic matrix

########################################    Models   #### ######################################################


#--------------------------#
# 2. Random Effects Models #
#--------------------------#


# Simple model (no random effects)
meta1 <- rma.uni(es, v, data= rdata, method= "REML")
summary(meta1) 

# Adding four random factors
meta2 <- rma.mv(es, v, random= list(~ 1|Species.phylo, ~ 1|species, ~ 1|study_code, ~1|obs), 
                R= list(Species.phylo = phylo_matrix), data= rdata, method= "REML")
summary(meta2)
i2_ml(meta2, method=c("ratio")) # Heterogeneity at each random factor level


# Accounting for non-independence of data points from the same experiment
# Assumes a correlation of 0.5 between effect sizes from the same experiment 
rdata$shared_control <- factor(rdata$Effect.size.code)
vcv_shared <- impute_covariance_matrix(vi=rdata$v, cluster = rdata$shared_control, r=0.5)

# Add new variance matrix into the mixed-effects meta-analysis model
meta3 <- rma.mv(es, VCV_shared, random= list(~ 1|Species.phylo, ~ 1|species, ~ 1|study_code, ~ 1|shared_control, ~1|obs), 
                R= list(Species.phylo = phylo_matrix), data= rdata, method= "REML")
summary(meta3)
i2_ml(meta3, method=c("ratio")) # Heterogeneity at each random factor level


## without phylogeny or species
meta4 <- rma.mv(es, VCV_shared, random= list(~ 1|study_code, ~ 1|shared_control, ~1|obs), data= rdata, method= "REML")
summary(meta4)
i2_ml(meta4, method=c("ratio")) # Heterogeneity at each random factor level


## without phylogeny but with shared control 
meta5 <- rma.mv(es, VCV_shared, random= list(~ 1|species, ~ 1|study_code, ~ 1|shared_control, ~1|obs), data= rdata, method= "REML")
summary(meta5)
i2_ml(meta5, method=c("ratio")) # Heterogeneity at each random factor level


## without phylogeny, species or study_code 
meta7 <- rma.mv(es, VCV_shared, random= list(~ 1|shared_control, ~1|obs), data= rdata, method= "REML")
summary(meta7)
i2_ml(meta7, method=c("ratio")) # Heterogeneity at each random factor level


## without phylogeny, species or shared_control 
meta8 <- rma.mv(es, VCV_shared, random= list(~ 1|study_code, ~1|obs), data= rdata, method= "REML")
summary(meta8)
i2_ml(meta8, method=c("ratio")) # Heterogeneity at each random factor level

##### meta6 (no species) is the best model out of the above, according to AIC

#-----------------------#
# 3. Meta-regressions   #
#-----------------------#


# Single categorical factor added as a fixed effect
AICs <- c()# Record AIC of single categorical factor fixed effect models

# warm/cool
meta_trait_warm <- rma.mv(es, VCV_shared,  mod= ~warm.cool, random= list(~ 1|study_code, ~ 1|shared_control, ~1|obs), data= rdata, method= "REML")
summary(meta_trait_warm)
aic <- AIC(meta_trait_warm)
AICs <- rbind(AICs, aic)

# warm/cool -1
meta_trait_warm_nointer <- rma.mv(es, VCV_shared,  mod= ~warm.cool-1, random= list(~ 1|study_code, ~ 1|shared_control, ~1|obs), data= rdata, method= "REML")
summary(meta_trait_warm_nointer)
aic <- AIC(meta_trait_warm_nointer)
AICs <- rbind(AICs, aic)

# ref temp
meta_trait_ref <- rma.mv(es, VCV_shared,  mod= ~reftemp, random= list(~ 1|study_code, ~ 1|shared_control, ~1|obs), data= rdata, method= "REML")
summary(meta_trait_ref)
aic <- AIC(meta_trait_ref)
AICs <- rbind(AICs, aic)

# treat temp
meta_trait_treattemp <- rma.mv(es, VCV_shared,  mod= ~treattemp, random= list(~ 1|study_code, ~ 1|shared_control, ~1|obs), data= rdata, method= "REML")
summary(meta_trait_treattemp)
aic <- AIC(meta_trait_treattemp)
AICs <- rbind(AICs, aic)

# treat temp^2
meta_trait_treat2 <- rma.mv(es, VCV_shared,  mod= ~ poly(treattemp, degree=2, raw=TRUE), random= list(~ 1|study_code, ~ 1|shared_control, ~1|obs), data= rdata, method= "REML")
summary(meta_trait_treat2)
aic <- AIC(meta_trait_treat2)
AICs <- rbind(AICs, aic)

# treat temp^3
meta_trait_treat3 <- rma.mv(es, VCV_shared,  mod= ~ poly(treattemp, degree=3, raw=TRUE), random= list(~ 1|study_code, ~ 1|shared_control, ~1|obs), data= rdata, method= "REML")
summary(meta_trait_treat3)
aic <- AIC(meta_trait_treat3)
AICs <- rbind(AICs, aic)

# diff temp
meta_trait_diff <- rma.mv(es, VCV_shared,  mod= ~diff, random= list(~ 1|study_code, ~ 1|shared_control, ~1|obs), data= rdata, method= "REML")
summary(meta_trait_diff)
aic <- AIC(meta_trait_diff)
AICs <- rbind(AICs, aic)

# all fixed effects
meta_trait_all <- rma.mv(es, VCV_shared,  mod= ~warm.cool + diff + treattemp + reftemp, random= list(~ 1|study_code, ~ 1|shared_control, ~1|obs), data= rdata, method= "REML")
summary(meta_trait_all)
aic <- AIC(meta_trait_all)
AICs <- rbind(AICs, aic)

### binned treatment max temperatures.
rdata$bin.temp <- c(NA)

rdata$bin.temp[which(rdata$treattemp >= 40)] <- ">40" 
rdata$bin.temp[which(rdata$treattemp >= 35 & rdata$treattemp <40)] <- "35-40" 
rdata$bin.temp[which(rdata$treattemp >= 30 & rdata$treattemp <35)] <- "30-35" 
rdata$bin.temp[which(rdata$treattemp >= 25 & rdata$treattemp <30)] <- "25-30" 
rdata$bin.temp[which(rdata$treattemp >= 20 & rdata$treattemp <25)] <- "20-25"
rdata$bin.temp[which(rdata$treattemp >= 15 & rdata$treattemp <20)] <- "15-20"
rdata$bin.temp[which(rdata$treattemp <15)] <- "<15" 

rdata$bin.temp <- factor(rdata$bin.temp)

levels(rdata$bin.temp)
table(rdata$bin.temp)

# binned temps
meta_trait_bintemp <- rma.mv(es, VCV_shared,  mod= ~bin.temp-1,  random= list(~ 1|study_code, ~ 1|shared_control, ~1|obs), data= rdata, method= "REML")
summary(meta_trait_bintemp)

