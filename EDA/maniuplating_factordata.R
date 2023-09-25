effectdata <- read.csv("Survival project all pairwise.es.csv")
rdata <- effectdata

rdata$Sex[which(rdata$Sex == "Both sexes")] <- "Both"
rdata$Sex[which(rdata$Sex == "Both Sexes")] <- "Both"



rdata$Sex.exposed[which(rdata$Sex.exposed == "Both sexes")] <- "Both"
rdata$Sex.exposed[which(rdata$Sex.exposed == "Both Sexes")] <- "Both"
rdata$Sex.exposed[which(rdata$Sex.exposed == "both")] <- "Both"
rdata$Sex.exposed[which(rdata$Sex.exposed == "female")] <- "Female"
rdata$Sex.exposed[which(rdata$Sex.exposed == "Females only")] <- "Female"
rdata$Sex.exposed[which(rdata$Sex.exposed == "male")] <- "Male"
rdata$Sex.exposed[which(rdata$Sex.exposed == "Males only")] <- "Male"


rdata$Agricultural.importance[which(rdata$Agricultural.importance == "no")] <- "No"


rdata$Country.of.origin[which(rdata$Country.of.origin == "Brasil")] <- "Brazil"


rdata$Exposure.duration[which(rdata$Exposure.duration == "1 - 5 days")] <- "1 to 5 days"
rdata$Exposure.duration[which(rdata$Exposure.duration == "Less than 24 hours")] <- "< 24 hours"
rdata$Exposure.duration[which(rdata$Exposure.duration == "Less than 24 hrs")] <- "< 24 hours"
rdata$Exposure.duration[which(rdata$Exposure.duration == "more than 5 days")] <- "> 5 days"
rdata$Exposure.duration[which(rdata$Exposure.duration == "More than 5 days")] <- "> 5 days"
rdata$Exposure.duration[which(rdata$Exposure.duration == "> 5 days")] <- "More than 5 days"
rdata$Exposure.duration[which(rdata$Exposure.duration == "Yes")] <- "More than 5 days"



rdata$Life.stage.of.animal[which(rdata$Life.stage.of.animal == "adult")] <- "Adult"
rdata$Life.stage.of.animal[which(rdata$Life.stage.of.animal == "Adults")] <- "Adult"
rdata$Life.stage.of.animal[which(rdata$Life.stage.of.animal == "Adults only")] <- "Adult"
rdata$Life.stage.of.animal[which(rdata$Life.stage.of.animal == "Eggs")] <- "Egg"
rdata$Life.stage.of.animal[which(rdata$Life.stage.of.animal == "from egg")] <- "Egg"
rdata$Life.stage.of.animal[which(rdata$Life.stage.of.animal == "from juvenile")] <- "Juvenile"
rdata$Life.stage.of.animal[which(rdata$Life.stage.of.animal == "Juveniles")] <- "Juvenile"
rdata$Life.stage.of.animal[which(rdata$Life.stage.of.animal == "Larval")] <- "Larvae"
rdata$Life.stage.of.animal[which(rdata$Life.stage.of.animal == "mix")] <- "Mix"
rdata$Life.stage.of.animal[which(rdata$Life.stage.of.animal == "No")] <- "Mix"


write.csv(rdata, "Survival project all pairwise.es.csv")
