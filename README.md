# Thermal_Fertility_Meta-Analysis

This folder contains everything needed to perform three seperate meta-analysis on how:

  - Temperature affects reproduction     (FILE: reproduction_meta-analysis_excHUM251.r)
  - Temperaature affects longevity       (FILE: longevity_meta-analysis_excHUM251.r)
  - Tempertature affects survival        (FILE: survival_meta-analysis_excHUM251.r)

NB: We have removed the study HUM251 in this uploaded analysis, hence the _excHUM251 extensions to many named files. 

  The code for each analysis initially loads the datafile containing the effect sizes:
  
"Survival project all pairwise.es.csv" 
  
  Then a specific subset of the data is taken depending which analysis you are performing i.e. only reproduction effect sizes are taken in the "reproduction_meta-analysis_excHUM251.r" file.

  The datafile 

"Species_Classifications.csv"

is called to enure that all species classifications are accurate.

  Depending on which analysis you are performing the specific phylogentic tree is loaded. This will be one of the three trees
   - all_long_excHUM251_tree.nex
   - all_reproduction_excHUM251_tree.nex
   - all_surv_excHUM251_tree.nex

   The code explores many models using combinations of random and fixed effects. 

   Since the meta-analysis is quite large and takes around 15minutes to compute each of the many models, the Rdata file for each analysis is uploaded containing the results. These files are named:

   - reproduction_meta-analysis_excHUM251.RData
   - longevity_meta-analysis.RData

     The Rdata fil for the survival meta-analysis is not provided as the number of effect sizes is quite small and the models can be run quickly on the fly.  
