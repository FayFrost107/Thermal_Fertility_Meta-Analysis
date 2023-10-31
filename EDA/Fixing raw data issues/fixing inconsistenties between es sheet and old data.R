geary <- read.csv("Data/Gearys_test_data.csv")
geary <- geary[,-1]

effectdata <- read.csv("Data/Survival project all pairwise.es.csv")


effectdata$trial.code <- paste(effectdata$Effect.size.code, effectdata$treattemp, sep="_")
geary$trial.code <- paste(geary$Effect.size.code, geary$treattemp, sep="_")



not_in_geary <- setdiff(effectdata$trial.code, geary$trial.code)
not_in_effectdata <- setdiff(geary$trial.code, effectdata$trial.code)



##############################
new_es <- read.csv("Survival project all pairwise.es.csv")# << from old folders 
current_es <- read.csv("Data/Survival project all pairwise.es.csv") # << from github



current_es$trial.code <- paste(current_es$Effect.size.code, current_es$treattemp, sep="_")
new_es$trial.code <- paste(new_es$Effect.size.code, new_es$treattemp, sep="_")


not_in_current <- setdiff(new_es$trial.code, current_es$trial.code)
not_in_new <- setdiff(current_es$trial.code, new_es$trial.code)

unique(not_in_current)

newdata <- new_es[!(new_es$trial.code %in% not_in_current), ]
