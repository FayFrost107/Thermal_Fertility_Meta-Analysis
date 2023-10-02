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

rdata <- rbind(repdata_warm, repdata_cool)

 

### Read in effect size data
longdata_warm <- subset(effectdata, Trait.category == "Longevity" & warm.cool == "Warm" )
longdata_cool <- subset(effectdata, Trait.category == "Longevity" & warm.cool == "Cool" )

alllong <- rbind(longdata_warm, longdata_cool)


rdata$es_reproduction <- rdata$es
rdata$v_reproduction <- rdata$v
rdata <- subset(rdata, select = -es)
rdata <- subset(rdata, select = -v)
rdata <- subset(rdata, select = -(Trait.category))
rdata <- subset(rdata, select = -(Trait))
rdata <- subset(rdata, select = -Effect.size.code)


alllong$es_longevity <- alllong$es
alllong$v_longevity <- alllong$v
alllong <- subset(alllong, select = -es)
alllong <- subset(alllong, select = -v)
alllong <- subset(alllong, select = -(Trait.category))
alllong <- subset(alllong, select = -(Trait))
alllong <- subset(alllong, select = -Effect.size.code)

df <- merge(rdata, alllong, all=TRUE)


write.csv(df, here("Data", "Combined_effectsizes.csv"))


df_cleaned <- df[!is.na(df$es_reproduction), ]
df_cleaned <- df_cleaned[!is.na(df_cleaned$es_longevity), ]


write.csv(df_cleaned, "Data/cleaned_unique_combo.csv")



########################## Find matching study info rows for reproduction and longevity dataframes and combine effect sizes #################
##########################  code below redundant. Have done this with merge function above. #####################################


#allrows <- rbind(rdata[,1:36], alllong[,1:36])

#combo_data <- allrows

#combo_data$es.reproduction <- NA
#combo_data$v.reproduction <- NA
#combo_data$es.longevity <- NA
#combo_data$v.longevity <- NA


#for (i in 1:nrow(combo_data)) {
  
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
  
#}

#u_combo <- unique(combo_data)

### calculate number of rows which have entries in both reproduction and longevity




#--------------#
# Long to wide #
#--------------#

# Unclear what the final long data should be from this code. You should be able to use pivot_wider if it has been set up correctly. If this doesn't work, another option is to: 1) subset the longevity; 2) subset the reproduction; 3) merge the two subsets together using left_join. You will first need to filter out the rows that are not in both subsets.

wide_dat <- read.csv(here("Data", "Combined_effectsizes.csv"))

long_dat <- wide_dat %>% 
              pivot_longer(cols = c(es_reproduction,  es_longevtiy), names_to = "es")   %>% 
              data.frame()











