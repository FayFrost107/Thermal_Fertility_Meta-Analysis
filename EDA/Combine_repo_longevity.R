#--------------#
# 1. Packages  #
#--------------#
install.packages("pacman")
pacman::p_load(tidyverse, here)

#--------------#
# 1. Setup     #
#--------------#

### Read in effect size data
effectdata <- read.csv(here("Data", "Survival project all pairwise.es.csv"))
effectdata <-  subset(effectdata, Paper.code != "HUM251") 


wide_dat <- pivot_wider(effectdata, values_from = c(es, v), names_from = Trait.category)  %>% data.frame()

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



library(dplyr)
library(tidyr)
df <- merge(rdata, alllong, all=TRUE)
df <- subset(df, select = -(Trait.category))
df <- subset(df, select = -(Trait))
df <- df[,-1]
df <- subset(df, select = -Effect.size.code)


write.csv(df, "Data/Combined_effectsizes.csv")



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

















