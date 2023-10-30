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
effectdata <- read.csv("Data/Survival project all pairwise.es.csv")
survdata_warm <- subset(effectdata, Trait.category == "Survival" & warm.cool == "Warm" )
survdata_cool <- subset(effectdata, Trait.category == "Survival" & warm.cool == "Cool" )

allsurv <- rbind(survdata_warm, survdata_cool)

### select data for analysis
rdata <- allsurv

rdata <- subset(rdata, Paper.code != "HUM251")
rdata <- subset(rdata, Paper.code != "OSM205")

### Species names which need changing for phylogeny. 
rdata$Species.latin[which(rdata$Species.latin == "Cosmocomoidea ashmeadi")]              <- "Gonatocerus ashmeadi"	
rdata$Species.latin[which(rdata$Species.latin == "Cosmocomoidea triguttata")]   <- "Gonatocerus triguttatus"
rdata$Species.latin[which(rdata$Species.latin == "Mythimna roseilinea")]            <-  "Mythimna albipuncta"
rdata$Species.latin[which(rdata$Species.latin == "Daphnia australis")]            <-  "Daphniopsis australis"

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
tree1 <- read.nexus("Phylogeny/all_surv_excHUM251_tree.nex")
tree_grafen = compute.brlen(tree1, method="Grafen", power=1)
tree <- ape::multi2di(tree_grafen, random = TRUE) 
phylo_matrix <- vcv(tree_grafen, cor=TRUE, model="Brownian") # Make phylogenetic matrix

# use a randomization approach to deal with polytomies. 
# Could this this approach or another detailed here: https://search.r-project.org/CRAN/refmans/RRphylo/html/fix.poly.html



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
VCV_shared <- matrix(0, nrow = dim(rdata)[1], ncol = dim(rdata)[1])
rownames(VCV_shared) <- rdata[, "es"]
colnames(VCV_shared) <- rdata[, "es"]
shared_coord <- which(rdata[, "shared_control"] %in% rdata[duplicated(rdata[, "shared_control"]), "shared_control"] == TRUE) 

#new_vcv <- impute_covariance_matrix(vi=rdata$v, cluster = rdata$study_code, r=0.5)


# Finds effect sizes that share a control group
combinations <- do.call("rbind", tapply(shared_coord, rdata[shared_coord,  "shared_control"], function(x) t(utils::combn(x, 2)))) 
for (i in 1:dim(combinations)[1]) {
  p1 <- combinations[i, 1]
  p2 <- combinations[i, 2]
  p1_p2_cov <- 0.5 * sqrt(rdata[p1, "v"]) * sqrt(rdata[p2, "v"])
  VCV_shared[p1, p2] <- p1_p2_cov
  VCV_shared[p2, p1] <- p1_p2_cov
} # Calculates the covariance between effect sizes and enters them in each combination of coordinates
diag(VCV_shared) <- rdata[, "v"] # Enters recalculated effect size sampling variances into diagonals


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

