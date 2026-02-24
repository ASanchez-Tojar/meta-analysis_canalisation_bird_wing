[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18152048.svg)](https://doi.org/10.5281/zenodo.18152048)

# Short title: Bird wing canalisation
## Title: Wing length canalisation and behaviour across birds: a phylogenetic meta-analysis of variance
**Contributors:** Klaus Reinhold, Alfredo Sánchez-Tójar 

**Date created:** March 2025

**Identifer:** [10.5281/zenodo.18152048](http://doi.org/10.5281/zenodo.18152048)

**Category:** Project

**Code License:** [MIT](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/blob/main/LICENSE)

**Data License:** [CC BY 4.0](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/blob/main/data/data_LICENSE.txt)

**Description:** This repository contains the data, code and other materials used in the following study:

---

Klaus Reinhold, Alfredo Sánchez-Tójar. 2026. Wing length canalisation and behaviour across birds: a phylogenetic meta-analysis of variance. Evolution (accepted). Preprint available at: https://doi.org/10.32942/X2PG92

---

The repository consists of an Rproject with 6 scripts and 4 folders containing the data either needed to run these scripts or created by these scripts, including the figures. For a detailed description of all the files, please see below. To run this repository, first download and unzip the entire repository and open the Rproject file:  body_size_variance_birds.Rproj. The scripts are named in the order that they should be run, from 001 to 005. For any further information about this repository, please contact: Alfredo Sánchez-Tójar, email: alfredo.tojar@gmail.com.

Information about folders and files within:
Folders:
*	[code](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/tree/main/code): contains all six scripts used to process and analyse data, as well as to run sensitivity analyses and generate the figures for the project.
    - [001_data_cleaning](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/blob/main/code/001_data_cleaning.R): used to process and clean our data.
    - [001B_The_new_phylogeny_after_review](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/blob/main/code/001B_The_new_phylogeny_after_review.R): this script was added after the first round of peer-review at the journal Evolution, during which we were asked to used a different, more updated bird tree. This script does that as well as updating the dataset with the new species column.
    - [002_data_analysis](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/blob/main/code/002_data_analysis.R): used to run the main analyses (as well as some of the sensitivity analyses).
    - [003_re-analysis_of_Nakagawa_et_al_2015](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/blob/main/code/003_re-analysis_of_Nakagawa_et_al_2015.R): used to re-analyse the data from Nakagawa et al. (2015).
    - [004_additional_sensitivity_analyses](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/blob/main/code/004_additional_sensitivity_analyses.R): used to run additional sensitivity analyses.
    - [005_phylogeny_figures](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/blob/main/code/005_phylogeny_figures.R): used to generate Figure 1 (phylogenetic tree with effect sizes).

*	[data](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/tree/main/data): contains all [final](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/tree/main/data/final) datasets used in the analyses as well as the [phylogeny](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/tree/main/data/phylogeny) data generated and used for the phylogenetic multilevel meta-analyses. The folder [Nakagawa_et_al_2015](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/tree/main/data/Nakagawa_et_al_2015) contains the data corresponding to the re-analysis of the data provided by Nakagawa et al. (2015).

*	[figures](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/tree/main/figures): contains all figures exported (including the supplementary figures).

*	[models](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/tree/main/models): contains all models run to faciliate reproducibility.

## Metadata:
We provide the metadata for the raw dataset [here](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/blob/main/data/final/METADATA.csv). Importantly, any new variable or file generated after the first round of peer-review at the journal Evolution is named with "_new" at the end.

## Software and {Packages} used
R version 4.3.1 (2023-06-16 ucrt)

Platform: x86_64-w64-mingw32/x64 (64-bit)

Running under: Windows 10 x64 (build 19044)

Packages attached: 
- rsvg_2.6.1          
- rgbif_3.8.1         
- ggcorrplot_0.1.4.1
- phytools_1.9-16
- maps_3.4.1.1
- readxl_1.4.3
- clootl_0.1.2
- clubSandwich_0.5.10
- orchaRd_2.0
- ggpubr_0.6.0       
- RColorBrewer_1.1-3  
- ggtree_3.10.1       
- ggstance_0.3.6      
- patchwork_1.3.0     
- viridis_0.6.5       
- viridisLite_0.4.2  
- cowplot_1.1.1       
- visdat_0.6.0        
- ggsignif_0.6.4      
- wesanderson_0.3.7   
- diagram_1.6.5       
- shape_1.4.6.1      
- treebase_0.1.4      
- metafor_4.6-0       
- numDeriv_2016.8-1.1
- metadat_1.2-0
- Matrix_1.6-1        
- ggstatsplot_0.12.5
- lubridate_1.9.4
- forcats_1.0.0       
- stringr_1.5.1       
- dplyr_1.1.4         
- purrr_1.0.2         
- readr_2.1.5        
- tidyr_1.3.1         
- tibble_3.2.1        
- ggplot2_3.5.1       
- tidyverse_2.0.0
- rotl_3.1.0
- ape_5.7-1

loaded via a namespace (and not attached):

  [1] mathjaxr_1.6-0          rstudioapi_0.17.1       jsonlite_1.8.8          datawizard_1.2.0       
  [5] correlation_0.8.6       magrittr_2.0.3          TH.data_1.1-3           estimability_1.5.1     
  [9] farver_2.1.2            fs_1.6.3                vctrs_0.6.5             memoise_2.0.1          
 [13] paletteer_1.6.0         RCurl_1.98-1.12         effectsize_0.8.9        rstatix_0.7.2          
 [17] progress_1.2.3          plotrix_3.8-2           broom_1.0.7             cellranger_1.1.0       
 [21] gridGraphics_0.5-1      plyr_1.8.9              sandwich_3.1-0          emmeans_1.10.6         
 [25] zoo_1.8-14              cachem_1.0.8            whisker_0.4.1           igraph_2.1.4           
 [29] lifecycle_1.0.4         iterators_1.0.14        pkgconfig_2.0.3         R6_2.6.1               
 [33] fastmap_1.1.1           digest_0.6.35           aplot_0.2.2             colorspace_2.1-0       
 [37] rematch2_2.1.2          clusterGeneration_1.3.8 timechange_0.3.0        httr_1.4.7             
 [41] abind_1.4-5             compiler_4.3.1          withr_3.0.2             doParallel_1.0.17      
 [45] backports_1.5.0         optimParallel_1.0-2     carData_3.0-5           MASS_7.3-60            
 [49] scatterplot3d_0.3-44    oai_0.4.0               tools_4.3.1             rncl_0.8.7             
 [53] rentrez_1.2.3           statsExpressions_1.6.1  glue_1.8.0              quadprog_1.5-8         
 [57] nlme_3.1-162            grid_4.3.1              generics_0.1.3          gtable_0.3.6           
 [61] tzdb_0.4.0              data.table_1.16.4       hms_1.1.3               xml2_1.4.0             
 [65] car_3.1-2               foreach_1.5.2           pillar_1.10.1           yulab.utils_0.1.4      
 [69] splines_4.3.1           treeio_1.26.0           lattice_0.21-8          survival_3.5-5         
 [73] tidyselect_1.2.1        gridExtra_2.3           expm_0.999-7            stringi_1.8.3          
 [77] lazyeval_0.2.2          ggfun_0.1.4             pacman_0.5.1            codetools_0.2-19       
 [81] ggplotify_0.1.2         cli_3.6.1               xtable_1.8-4            parameters_0.23.0      
 [85] munsell_0.5.1           Rcpp_1.1.0              zeallot_0.1.0           coda_0.19-4            
 [89] XML_3.99-0.14           parallel_4.3.1          prettyunits_1.2.0       bayestestR_0.15.0      
 [93] bitops_1.0-7            phangorn_2.12.1         mvtnorm_1.3-2           tidytree_0.4.6         
 [97] scales_1.3.0            insight_1.4.2           crayon_1.5.3            combinat_0.0-8         
[101] rlang_1.1.4             fastmatch_1.1-4         multcomp_1.4-28         mnormt_2.1.1 
