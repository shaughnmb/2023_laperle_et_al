# 2022_laperle_et_al
R script for the publication:  "iPSC-derived neural progenitor cells secreting GDNF slow disease progression in rodent models of both ALS and Retinal degeneration"


Title: iPSC-derived neural progenitor cells secreting GDNF slow disease
        progression in rodent models of both ALS and Retinal degeneration
Authors:  Alexander H. Laperle* 1, Alexandra Moser* 1, Veronica J. Garcia 1,
          Amanda Wu 1, Aaron Fulton 1, George Lawless 1, Shaughn Bell 1,
          Kristina Roxas 1, Roksana Elder 1, Pablo Avalos 1, Bin Lu 1,
          Staphany Ramirez 1, Shaomei Wang 1, Clive N. Svendsen** 1

Affiliations: 1 Cedars-Sinai Board of Governors Regenerative Medicine
              Institute, Cedars-Sinai Medical Center, Los Angeles CA
 '* These authors contributed equally to this work
 ** Corresponding author email:  Clive.Svendsen@cshs.org


R version 4.2.0
R Script Title:  laperele_et_al_2022_v2.R
R Script Author:  Shaughn Bell
R Script Corresponding Email:  shaughn.bell@cshs.org

Notes: 
   A) Script makes use of the variables set up under "project information" as 
      well as additional "prefixes" throughout the script for ease of saving
      files with a similar path and naming structure.  When reloading data 
      (see note "B"), you must either use the full path or reload the prefixes.
   B) Script saves intermediate steps at each major manipulation of the seurat
      object via the "saveRDS" function.  If needed, these RDS objects can then 
      be reloaded to easily restart at one of these save points without needing 
      to start from scratch.  However, these are not required for analysis, and
      they can be skipped to save time and disk space.
