library(tidyr)
library(dplyr)

# get effectsize data
effectdata <- read.csv("Data/Survival project all pairwise.es.csv")
sub <- subset(effectdata, warm.cool != "Reference")

# How many papers in total
nlevels(as.factor(sub$Paper.code))

## percentage effect sizes of each class 
table(sub$Class)/nrow(sub) * 100


insect <- subset(sub, Class == "Insecta")
arach <- subset(sub, Class == "Arachnida")
crus <- subset(sub, Class == "Crustacea")
rot <- subset(sub, Class == "Rotifera")
ann <- subset(sub, Class == "Annelida")

nlevels(as.factor(insect$Species.latin))
nlevels(as.factor(arach$Species.latin))
nlevels(as.factor(crus$Species.latin))
nlevels(as.factor(rot$Species.latin))
nlevels(as.factor(ann$Species.latin))

# How many papers does each species occur in
unique_combinations <- sub %>%
  distinct(Paper.code, Species.latin)

# Count the number of unique studies for each species
species_counts <- unique_combinations %>%
  group_by(Species.latin) %>%
  summarize(UniqueStudies = n())

# number of species
length(unique(unique_combinations$Species.latin))


# Terrestiral versus aquatic
unique_habitat <- sub %>%
   distinct(Species.latin, Habitat)
 
# Count the number of unique studies for each species
habitat_counts <- unique_habitat %>%
     group_by(Habitat) %>%
     summarize(UniqueStudies = n())

habitat_counts

unique_habitat2 <- sub %>%
  distinct(Species.latin, Habitat2)

# Count the number of unique studies for each species
habitat2_counts <- unique_habitat2 %>%
  group_by(Habitat2) %>%
  summarize(UniqueStudies = n())

habitat2_counts

### Fertiliarsation mode

# internal/external
unique_fert <- sub %>%
  distinct(Species.latin, Fertilisation.mode)

# Count the number of unique studies for each species
fert_counts <- unique_fert %>%
  group_by(Fertilisation.mode) %>%
  summarize(UniqueStudies = n())

fert_counts


### reproduction mode

# sex/asexual
unique_sex <- sub %>%
  distinct(Species.latin, reprodctuive.mode)

# Count the number of unique studies for each species
sex_counts <- unique_sex %>%
  group_by(reprodctuive.mode) %>%
  summarize(UniqueStudies = n())

sex_counts

### lab per studies

# lab study
unique_lab <- sub %>%
  distinct(Experiment.code, Lab.or.field)

# Count the number of lab studies for each study
lab_counts <- unique_lab %>%
  group_by(Lab.or.field) %>%
  summarize(UniqueStudies = n())

lab_counts
# exposure
unique_exp <- sub %>%
  distinct(Experiment.code, Exposure.duration)

# Count the number of exposure times for each study
exp_counts <- unique_exp %>%
  group_by(Exposure.duration) %>%
  summarize(UniqueStudies = n())

exp_counts


