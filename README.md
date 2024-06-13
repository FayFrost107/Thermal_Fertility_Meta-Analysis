# Thermal_Fertility_Meta-Analysis

This folder contains everything needed to perform three seperate meta-analysis on how:

  - Temperature affects reproduction     (R/ meta_analysis_reproduction.r)
  - Temperaature affects longevity       (R/ meta_analysis_longevity.r)
  - Tempertature affects survival        (R/ meta_analysis_survival.r)

NB: We have removed the study HUM251 in this uploaded analysis, hence the _excHUM251 extensions to many named files. 

  The code for each analysis initially loads the datafile containing the effect sizes:
  
Data / "Survival project all pairwise.es.csv" 
  
  Then a specific subset of the data is taken depending which analysis you are performing i.e. only reproduction effect sizes are taken in the "reproduction_meta-analysis_excHUM251.r" file.

  The datafile 

Data / "Species_Classifications.csv"

is called to enure that all species classifications are accurate.


 Depending on which analysis you are performing the specific phylogentic tree is loaded. This will be one of the three trees in the phylogeny folder
  
   - all_long_excHUM251_tree.nex
   - all_reproduction_excHUM251_tree.nex
   - all_surv_excHUM251_tree.nex

   The code explores many models using combinations of random and fixed effects. 

   Since the meta-analysis is quite large and takes around 15minutes to compute each of the many models, the Rdata file for each analysis is uploaded containing the results. These files are named:

   - Output / reproduction_meta-analysis_excHUM251.RData
   - Output / longevity_meta-analysis.RData

     The Rdata file for the survival meta-analysis is not provided as the number of effect sizes is quite small and the models can be run quickly on the fly.

A summary of the completed analysis for temperature effects on reproduction and longevity are given in the Rmarkdown folder. The files "meta_reproduction.qmd" and "meta_longevity.qmd" provide a comprehensive explanation of the final models used in each analysis and the results that were obtained. 
