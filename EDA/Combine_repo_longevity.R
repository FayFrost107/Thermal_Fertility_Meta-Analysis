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



#################### trying to combne repro and longevity below #### currently not working

### Read in effect size data
effectdata <- read.csv("Data/Survival project all pairwise.es.csv")
effectdata <-  subset(effectdata, Paper.code != "HUM251") 


repdata_warm <- subset(effectdata, Trait.category == "Reproduction" & warm.cool == "Warm" )
repdata_cool <- subset(effectdata, Trait.category == "Reproduction" & warm.cool == "Cool" )

allrep <- rbind(repdata_warm, repdata_cool)

rdata <- allrep
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

alllong <- subset(alllong, select = -(Trait.category))
alllong <- subset(alllong, select = -(Trait))
alllong <- alllong[,-1]
alllong <- subset(alllong, select = -Effect.size.code)

rdata <- subset(rdata, select = -(Trait.category))
rdata <- subset(rdata, select = -(Trait))
rdata <- rdata[,-1]
rdata <- subset(rdata, select = -Effect.size.code)

allrows <- rbind(rdata[,1:35], alllong[,1:35])
unique_rows <- allrows[!duplicated(allrows), ]

combo_data <- unique_rows



for(i in 1:nrow(unique_rows)){

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

 
    combo_data$es.reproduction[i] <- rdata$es_reproduction[is.identical_r]
    combo_data$v.reproduction[i] <-  rdata$v_reproduction[is.identical_r]
    
    combo_data$es.longevity[i] <- alllong$es_longevity[is.identical_l]
    combo_data$v.longevity[i] <-  alllong$v_longevity[is.identical_l]

}



#--------------#
# Long to wide #
#--------------#

# Unclear what the final long data should be from this code. You should be able to use pivot_wider if it has been set up correctly. If this doesn't work, another option is to: 1) subset the longevity; 2) subset the reproduction; 3) merge the two subsets together using left_join. You will first need to filter out the rows that are not in both subsets.

# Note that I am not sure what the final long_data file is, is it the "combined_effectsizes.csv"? Assuming it's called "long_data" you can use the code below to pivot_wider.
wide_dat <- long_data %>% 
              pivot_wider(values_from = c(es, v), names_from = Trait.category)   %>% 
              data.frame()











