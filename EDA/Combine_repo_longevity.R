

#--------------#
# 1. Setup     #
#--------------#

### Read in effect size data
effectdata <- read.csv("Data/Survival project all pairwise.es.csv")
repdata_warm <- subset(effectdata, Trait.category == "Reproduction" & warm.cool == "Warm" )
repdata_cool <- subset(effectdata, Trait.category == "Reproduction" & warm.cool == "Cool" )

allrep <- rbind(repdata_warm, repdata_cool)

rdata <- allrep

rdata <- subset(rdata, Paper.code != "HUM251")

### Read in effect size data
effectdata <- read.csv("Data/Survival project all pairwise.es.csv")
longdata_warm <- subset(effectdata, Trait.category == "Longevity" & warm.cool == "Warm" )
longdata_cool <- subset(effectdata, Trait.category == "Longevity" & warm.cool == "Cool" )

alllong <- rbind(longdata_warm, longdata_cool)

alllong <- subset(alllong, Paper.code != "HUM251")

rdata$es_reproduction <- rdata$es
rdata$v_reproduction <- rdata$v

rdata <- subset(rdata, select = -es)
rdata <- subset(rdata, select = -v)


alllong$es_longevtiy <- alllong$es
alllong$v_longevtiy <- alllong$v
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
