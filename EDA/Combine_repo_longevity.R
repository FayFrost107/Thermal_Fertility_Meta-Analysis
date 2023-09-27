#--------------#
# 1. Packages  #
#--------------#
install.packages("pacman")
pacman::p_load(tidyverse, here, dplyr, tidyr)


#--------------#
# 1. Setup     #
#--------------#

### Read in effect size data
effectdata <- read.csv(here("Data", "Survival project all pairwise.es.csv"))
effectdata <-  subset(effectdata, Paper.code != "HUM251") 


repdata_warm <- subset(effectdata, Trait.category == "Reproduction" & warm.cool == "Warm" )
repdata_cool <- subset(effectdata, Trait.category == "Reproduction" & warm.cool == "Cool" )

allrep <- rbind(repdata_warm, repdata_cool)


### Read in effect size data
longdata_warm <- subset(effectdata, Trait.category == "Longevity" & warm.cool == "Warm" )
longdata_cool <- subset(effectdata, Trait.category == "Longevity" & warm.cool == "Cool" )

alllong <- rbind(longdata_warm, longdata_cool)


rdata$es_reproduction <- rdata$es
rdata$v_reproduction <- rdata$v

rdata <- subset(rdata, select = -es)
rdata <- subset(rdata, select = -v)


alllong$es_longevity <- alllong$es
alllong$v_longevity <- alllong$v
alllong <- subset(alllong, select = -es)
alllong <- subset(alllong, select = -v)


df <- merge(rdata, alllong, all=TRUE)
df <- subset(df, select = -(Trait.category))
df <- subset(df, select = -(Trait))
df <- df[,-1]
df <- subset(df, select = -Effect.size.code)


write.csv(df, here("Data", "Combined_effectsizes.csv"))



#################### trying to combine repro and longevity below #### currently not working

### Read in effect size data
effectdata <- read.csv("Data/Survival project all pairwise.es.csv")
effectdata <-  subset(effectdata, Paper.code != "HUM251") ## remove outlier paper

## subset reproduction data
repdata_warm <- subset(effectdata, Trait.category == "Reproduction" & warm.cool == "Warm" )
repdata_cool <- subset(effectdata, Trait.category == "Reproduction" & warm.cool == "Cool" )

allrep <- rbind(repdata_warm, repdata_cool)

rdata <- allrep


### subset longevity data
longdata_warm <- subset(effectdata, Trait.category == "Longevity" & warm.cool == "Warm" )
longdata_cool <- subset(effectdata, Trait.category == "Longevity" & warm.cool == "Cool" )

alllong <- rbind(longdata_warm, longdata_cool)

# make new variables
rdata$es_reproduction <- rdata$es
rdata$v_reproduction <- rdata$v

# remove old variables
rdata <- subset(rdata, select = -es)
rdata <- subset(rdata, select = -v)

# make new longevity variables and remove old
alllong$es_longevity <- alllong$es
alllong$v_longevity <- alllong$v
alllong <- subset(alllong, select = -es)
alllong <- subset(alllong, select = -v)


# remove columns which are now redundant
alllong <- subset(alllong, select = -(Trait.category))
alllong <- subset(alllong, select = -(Trait))
alllong <- alllong[,-1]
alllong <- subset(alllong, select = -Effect.size.code)

rdata <- subset(rdata, select = -(Trait.category))
rdata <- subset(rdata, select = -(Trait))
rdata <- rdata[,-1]
rdata <- subset(rdata, select = -Effect.size.code)

## put reproduction and longevity together
allrows <- rbind(rdata[,1:35], alllong[,1:35])

# find unique rows
unique_rows <- allrows[!duplicated(allrows), ]


combo_data <- unique_rows

combo_data$es.reproduction <- NA
combo_data$v.reproduction <- NA
combo_data$es.longevity <- NA
combo_data$v.longevity <- NA

for(i in 1:nrow(unique_rows)){
  
  cat("Calculating row", i, "/2062", "\n")
  
  is.identical_r <- c()
  is.identical_l <- c()
  
  for(j in 1:nrow(rdata)){
    
    result_r <- identical(rdata[j, 1:35], unique_rows[i,])
    is.identical_r <- c(result_r, is.identical_r)
  }
  
  for(j in 1:nrow(alllong)){
    result_l <- identical(alllong[j, 1:35], unique_rows[i,])
    is.identical_l <- c(result_l, is.identical_l)
  }
  
  if(any(is.identical_r == TRUE)){  ## if there is an true index then there are rep values
    combo_data$es.reproduction[i] <- rdata$es_reproduction[is.identical_r]
    combo_data$v.reproduction[i] <-  rdata$v_reproduction[is.identical_r]
  } 
  
  if(any(is.identical_l == TRUE)){ ## if there is an true index then there are long values
    combo_data$es.longevity[i] <- alllong$es_longevity[is.identical_l]
    combo_data$v.longevity[i] <-  alllong$v_longevity[is.identical_l]
  } 
  
  }


######################################################
## I think the above is doing the wrong thing as I am just looking at the unique rows there. 
## In other words, rows that just appear once in the data. 
## These won't have a longevity AND rep component, they'll just have one. 
## I need to keep exactly one of EVERY row, duplicated or not.  
## Maybe i should do it for allrows then delete any duplicates at the end .....
 

#### i think there are only 742 rows for which there is both longevity and repro

allrows <- rbind(rdata[,1:35], alllong[,1:35])

combo_data <- allrows

combo_data$es.reproduction <- NA
combo_data$v.reproduction <- NA
combo_data$es.longevity <- NA
combo_data$v.longevity <- NA


for (i in 1:nrow(combo_data)) {
  
  cat("Calculating row", i, "/2804", "\n")
  
  is_matching_r <- rep(FALSE, nrow(rdata))
  is_matching_l <- rep(FALSE, nrow(alllong))
  
  # for every row of allrows find matching row in rdata and alllong
  
  for (j in 1:nrow(rdata)) {
    if (all(allrows[i, 1:35] == rdata[j, 1:35], na.rm=T)) {
      is_matching_r[j] <- TRUE
      break  # Break out of the loop once a match is found
    }
  }
  
  for (j in 1:nrow(alllong)) {
    if (all(allrows[i, 1:35] == alllong[j, 1:35], na.rm=T)) {
      is_matching_l[j] <- TRUE
      break  # Break out of the loop once a match is found
    }
  }
  
  if (any(is_matching_r)) {  
    # Extract values from rdata and assign to combo_data
    combo_data$es.reproduction[i] <- rdata$es_reproduction[is_matching_r]
    combo_data$v.reproduction[i] <- rdata$v_reproduction[is_matching_r]
  } 
  
  if (any(is_matching_l)) { 
    # Extract values from alllong and assign to combo_data
    combo_data$es.longevity[i] <- alllong$es_longevity[is_matching_l]
    combo_data$v.longevity[i] <- alllong$v_longevity[is_matching_l]
  } 
  
}

#--------------#
# Long to wide #
#--------------#

# Unclear what the final long data should be from this code. You should be able to use pivot_wider if it has been set up correctly. If this doesn't work, another option is to: 1) subset the longevity; 2) subset the reproduction; 3) merge the two subsets together using left_join. You will first need to filter out the rows that are not in both subsets.

wide_dat <- read.csv(here("Data", "Combined_effectsizes.csv"))

long_dat <- wide_dat %>% 
              pivot_longer(cols = c(es_reproduction,  es_longevtiy), names_to = "es")   %>% 
              data.frame()











