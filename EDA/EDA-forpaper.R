# get effectsize data
effectdata <- read.csv("Data/Survival project all pairwise.es.csv")
sub <- subset(effectdata, warm.cool != "Reference")


########### change species names in survival data ####################################
classes <- read.csv("Data/Species_classifications.CSV") ## read in species classifications from map

sub$Species.latin[which(sub$Species.latin == "Marasmia exigua")]                <- "Cnaphalocrocis exigua"
sub$Species.latin[which(sub$Species.latin == "Matsumuratettix hieroglyphicus")] <- "Matsumuratettix hiroglyphicus"
sub$Species.latin[which(sub$Species.latin == "Mythimna roseilinea")]            <- "Mythimna albipuncta"
sub$Species.latin[which(sub$Species.latin == "Apis craccivora")]                <- "Aphis craccivora"
sub$Species.latin[which(sub$Species.latin == "Cryptoleamus montrouzieri")]      <- "Cryptolaemus montrouzieri"
sub$Species.latin[which(sub$Species.latin == "Asplanchna brightwelli")]         <- "Asplanchna brightwellii"
sub$Species.latin[which(sub$Species.latin == "Brennandania lambi")]             <- "Pygmephorus lambi"
sub$Species.latin[which(sub$Species.latin == "Amblyseius alstoniae")]           <- "Euseius alstoniae"
sub$Species.latin[which(sub$Species.latin == "Siphoninus phyllyreae")]          <- "Siphoninus phillyreae"
sub$Species.latin[which(sub$Species.latin == "Proprioseiopsis asetus")]         <- "Amblyseius asetus"
sub$Species.latin[which(sub$Species.latin == "Parabemisia myrica")]             <- "Parabemisia myricae"
sub$Species.latin[which(sub$Species.latin == "Cirrospilus sp. near lyncus")]    <- "Cirrospilus lyncus"
sub$Species.latin[which(sub$Species.latin == "Anagyrus sp. nov. nr. sinope" )]  <- "Anagyrus sinope"
sub$Species.latin[which(sub$Species.latin == "Monochamus leuconotus")]          <- "Anthores leuconotus"
sub$Species.latin[which(sub$Species.latin == "Ropalosiphum maidis")]            <- "Rhopalosiphum maidis"

sub$Species.latin[which(sub$Species.latin == "Artemia fransiscana")]            <- "Artemia franciscana"
sub$Species.latin[which(sub$Species.latin == "Blathyplectes curculionis")]      <- "Bathyplectes curculionis"
sub$Species.latin[which(sub$Species.latin == "Menochilus sexmaculatus")]        <- "Cheilomenes sexmaculata"
sub$Species.latin[which(sub$Species.latin == "unknown (Tominic)")]              <- "Trichogramma" 
### specify classifications from map
sub$Class <- classes$class[match(sub$Species.latin, classes$species_latin)]


## number of each class 
table(sub$Class)

insect <- subset(sub, Class == "Insecta")
arach <- subset(sub, Class == "Arachnida")
crus <- subset(sub, Class == "Crustacea")
rot <- subset(sub, Class == "Rotifera")
ann <- subset(sub, Class == "Annelida")


# How many papers does each species occur in
unique_combinations <- sub %>%
  distinct(Paper.code, Species.latin)

# Count the number of unique studies for each species
species_counts <- unique_combinations %>%
  group_by(Species.latin) %>%
  summarize(UniqueStudies = n())


# Terrestiral versus aquatic
unique_habitat <- sub %>%
   distinct(Species.latin, Habitat)
 
# Count the number of unique studies for each species
habitat_counts <- unique_habitat %>%
     group_by(Habitat) %>%
     summarize(UniqueStudies = n())


### Fertiliarsation mode

# internal/external
unique_fert <- sub %>%
  distinct(Species.latin, Fertilisation.mode)

# Count the number of unique studies for each species
fert_counts <- unique_fert %>%
  group_by(Fertilisation.mode) %>%
  summarize(UniqueStudies = n())



### lab per studies

# lab study
unique_lab <- sub %>%
  distinct(Experiment.code, Lab.or.field)

# Count the number of lab studies for each study
lab_counts <- unique_lab %>%
  group_by(Lab.or.field) %>%
  summarize(UniqueStudies = n())


# exposure
unique_exp <- sub %>%
  distinct(Experiment.code, Exposure.duration)

# Count the number of exposure times for each study
exp_counts <- unique_exp %>%
  group_by(Exposure.duration) %>%
  summarize(UniqueStudies = n())




