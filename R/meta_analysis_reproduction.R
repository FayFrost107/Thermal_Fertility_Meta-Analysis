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
library(tidyr)

# To install the orchaRd package:
#install.packages("pacman")
#pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, emmeans)
#devtools::install_github("daniel1noble/orchaRd", force = TRUE)
library(orchaRd)
library(ggtree)
library(clubSandwich)
library(rmarkdown)

##########################################################################################################


#--------------#
# 1. Setup     #
#--------------#

### Read in effect size data
effectdata <- read.csv("Data/Survival project all pairwise.es.csv")

### select data for analysis
repdata_warm <- subset(effectdata, Trait.category == "Reproduction" & warm.cool == "Warm" )
repdata_cool <- subset(effectdata, Trait.category == "Reproduction" & warm.cool == "Cool" )
allrep <- rbind(repdata_warm, repdata_cool)
rdata <- allrep

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
tree1 <- read.nexus("Phylogeny/all_rep_excHUM251_tree.nex")
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
VCV_shared <- impute_covariance_matrix(vi=rdata$v, cluster = rdata$shared_control, r=0.5)


# Add new variance matrix for shared_control into the mixed-effects meta-analysis model
meta3 <- rma.mv(es, VCV_shared, random= list(~ 1|Species.phylo, ~ 1|species, ~ 1|study_code, ~1|obs), 
                R= list(Species.phylo = phylo_matrix), data= rdata, method= "REML")
summary(meta3)
i2_ml(meta3, method=c("ratio")) # Heterogeneity at each random factor level


## without phylogeny or species
meta4 <- rma.mv(es, VCV_shared, random= list(~ 1|study_code, ~1|obs), data= rdata, method= "REML")
summary(meta4)
i2_ml(meta4, method=c("ratio")) 


## without phylogeny 
meta5 <- rma.mv(es, VCV_shared, random= list(~ 1|species, ~ 1|study_code,  ~1|obs), data= rdata, method= "REML")
summary(meta5)
i2_ml(meta5, method=c("ratio")) 


## without phylogeny, species or study_code 
meta7 <- rma.mv(es, VCV_shared, random= list(~1|obs), data= rdata, method= "REML")
summary(meta7)
i2_ml(meta7, method=c("ratio")) # Heterogeneity at each random factor level


#-----------------------#
# 3. Meta-regressions   #
#-----------------------#


# Single categorical factor added as a fixed effect

# warm/cool
meta_trait_warm <- rma.mv(es, VCV_shared,  mod= ~warm.cool, random= list(~ 1|study_code,  ~1|obs), data= rdata, method= "REML")
summary(meta_trait_warm)

# warm/cool -1
meta_trait_warm_nointer <- rma.mv(es, VCV_shared,  mod= ~warm.cool-1, random= list(~ 1|study_code,  ~1|obs), data= rdata, method= "REML")
summary(meta_trait_warm_nointer)

# ref temp
meta_trait_ref <- rma.mv(es, VCV_shared,  mod= ~reftemp, random= list(~ 1|study_code,  ~1|obs), data= rdata, method= "REML")
summary(meta_trait_ref)


# treat temp centered
meta_trait_treattemp <- rma.mv(es, VCV_shared,  mod= ~c_treattemp, random= list(~ 1|study_code,  ~1|obs), data= rdata, method= "REML")
summary(meta_trait_treattemp)

# treat temp^2 centered
meta_trait_treat2 <- rma.mv(es, VCV_shared,  mod= ~poly(c_treattemp, degree=2, raw=TRUE), random= list(~ 1|study_code,  ~1|obs), data= rdata, method= "REML")
summary(meta_trait_treat2)

saveRDS(meta_trait_treat2, here("output", "models", "meta_rep_2.rds"))

# diff temp
meta_trait_diff <- rma.mv(es, VCV_shared,  mod= ~diff, random= list(~ 1|study_code,  ~1|obs), data= rdata, method= "REML")
summary(meta_trait_diff)


# all fixed effects
meta_trait_all <- rma.mv(es, VCV_shared,  mod= ~warm.cool + diff + treattemp + reftemp, random= list(~ 1|study_code,  ~1|obs), data= rdata, method= "REML")
summary(meta_trait_all)


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
meta_trait_bintemp <- rma.mv(es, VCV_shared,  mod= ~bin.temp-1,  
                             random= list(~ 1|study_code,  ~1|obs), data= rdata, method= "REML")
summary(meta_trait_bintemp)



## Publication Bias.

meta_year <- rma.mv(es, VCV_shared,  mod= ~Publication.year,  
                    random= list(~ 1|study_code,  ~1|obs), data= rdata, method= "REML")
summary(meta_year)


# Sensitivty Analysis
# Preform a sensitivity analysis by removing the smallest and largest 2.5% of effect sizes.


minq <- quantile(rdata$es, 0.025)
maxq <- quantile(rdata$es, 0.975)

sdata <- subset(rdata, es > minq & es < maxq)

## Treatment temperature as a cubic effect (sesnsitivity analysis)

# re-cmpute the covariance matrix for subsetted data
VCV_shared_sa <- impute_covariance_matrix(vi=sdata$v, cluster = sdata$shared_control, r=0.5)


# re-run quadratic model
meta_sa_treat2 <- rma.mv(es, VCV_shared_sa,  mod= ~poly(c_treattemp, degree=2, raw=TRUE), 
                         random= list(~ 1|study_code,  ~1|obs), data= sdata, method= "REML")
summary(meta_sa_treat2)


# re-run binned temps model
meta_sa_bintemp <- rma.mv(es, VCV_shared_sa,  mod= ~bin.temp-1,  
                          random= list(~ 1|study_code,  ~1|obs), data= sdata, method= "REML")

summary(meta_sa_bintemp)


## Now we completely remove any study that has an effect size in the highest or lowest 2.5%. 

remove_min <- unique(rdata$Paper.code[which(rdata$es < minq)])
remove_max <- unique(rdata$Paper.code[which(rdata$es > maxq)])
remove <- union(remove_min, remove_max)

subdata <- rdata[!(rdata$Paper.code %in% remove), ]

VCV_shared_subdata <- impute_covariance_matrix(vi=subdata$v, cluster = subdata$shared_control, r=0.5)

meta_sub_bintemp <- rma.mv(es, VCV_shared_subdata,  mod= ~bin.temp-1,  
                           random= list(~ 1|study_code,  ~1|obs), data= subdata, method= "REML")

summary(meta_sub_bintemp)

###########################################################################################################
#Other fixed effects


## Sex exposed
# We could lump categories so that we have cases where males are included (Both, Male), versus cases with just females (Female, Parthenogenetic), 
# with Unsure removed. I would predict that the 'Both' category would show the biggest drop for reproduction, but there will be no difference for lifespan

new_data <- rdata

new_data$Sex.exposed[which(new_data$Sex.exposed == "Male")] <- "Both"
new_data$Sex.exposed[which(new_data$Sex.exposed == "Parthenogenetic")] <- "Female"

new_data <- subset(new_data, Sex.exposed != "Unsure")

VCV_shared_sex <- impute_covariance_matrix(vi=new_data$v, cluster = new_data$shared_control, r=0.5)

meta_treat_sex <- rma.mv(es, VCV_shared_sex,  mod= ~poly(c_treattemp, degree=2, raw=TRUE)*Sex.exposed-1, 
                         random= list(~ 1|study_code,  ~1|obs), data= new_data, method= "REML")

saveRDS(meta_treat_sex, here("output", "models", "meta_treat_rep_sex.rds"))
saveRDS(new_data, here("output", "Output data", "data_rep_sex.rds"))

## Life-stage
# We could lump categories so that we have cases where only adults were exposed (Adult), 
# versus cases where immature stages were exposed (Juvenile, Larvae, Pupae, Mix)- perhaps after excluding 'Egg' and 'Embryo' because these categories are a bit weird




ls_data <- rdata

ls_data$Life.stage.of.animal[which(ls_data$Life.stage.of.animal == "Juvenile")] <- "Immature"
ls_data$Life.stage.of.animal[which(ls_data$Life.stage.of.animal == "Larvae")] <- "Immature"
ls_data$Life.stage.of.animal[which(ls_data$Life.stage.of.animal == "Mix")] <- "Immature"
ls_data$Life.stage.of.animal[which(ls_data$Life.stage.of.animal == "Pupae")] <- "Immature"



ls_data <- subset(ls_data, Life.stage.of.animal != "Egg")
ls_data <- subset(ls_data, Life.stage.of.animal != "Embryo")

VCV_shared_life <- impute_covariance_matrix(vi=ls_data$v, cluster = ls_data$shared_control, r=0.5)

meta_treat_ls <- rma.mv(es, VCV_shared_life,  mod= ~poly(c_treattemp, degree=2, raw=TRUE)*Life.stage.of.animal,
                        random= list(~ 1|study_code,  ~1|obs), data= ls_data, method= "REML")



saveRDS(meta_treat_ls, here("output", "models", "meta_treat_rep_ls.rds"))
saveRDS(ls_data, here("output", "Output data", "data_rep_ls.rds"))


### pest
pest_data <- subset(rdata, Agricultural.importance == "Pest")


VCV_shared_pest <- impute_covariance_matrix(vi=pest_data$v, cluster = pest_data$shared_control, r=0.5)

meta_pest <- rma.mv(es, VCV_shared_pest,  mod= ~poly(c_treattemp, degree=2, raw=TRUE),
                        random= list(~ 1|study_code,  ~1|obs), data= pest_data, method= "REML")
