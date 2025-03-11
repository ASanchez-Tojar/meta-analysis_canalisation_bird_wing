# Short title: Bird wing canalisation
## Title: Wing length canalisation and behaviour across birds: a phylogenetic meta-analysis of variance
**Contributors:** Klaus Reinhold, Alfredo Sánchez-Tójar 

**Date created:** February 2025

**Identifer:** TBA

**Category:** Project

**Code License:** [MIT](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/blob/main/LICENSE)

**Data License:** [CC BY 4.0](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/blob/main/data/data_LICENSE.txt)

**Description:** This repository contains the data, code and other materials used in the following study:

---

Klaus Reinhold, Alfredo Sánchez-Tójar. 2025. Wing length canalisation and behaviour across birds: a phylogenetic meta-analysis of variance. Preprint available at: TBA

---

The repository consists of an Rproject with 5 scripts and 4 folders containing the data either needed to run these scripts or created by these scripts, including the figures. For a detailed description of all the files, please see below. To run this repository, first download and unzip the entire repository and open the Rproject file:  body_size_variance_birds.Rproj. The scripts are named in the order that they should be run, from 001 to 005. For any further information about this repository, please contact: Alfredo Sánchez-Tójar, email: alfredo.tojar@gmail.com.

Information about folders and files within:
Folders:
*	[code](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/tree/main/code): contains all five scripts used to process and analyse data, as well as to run sensitivity analyses and generate the figures for the project.
    - [001_data_cleaning](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/blob/main/code/001_data_cleaning): used to process and clean our data.
    -	[002_data_analysis](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/blob/main/code/002_data_analysis): used to run the main analyses (as well as some of the sensitivity analyses).
    -	[003_re-analysis_of_Nakagawa_et_al_2015](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/blob/main/code/003_re-analysis_of_Nakagawa_et_al_2015): used to re-analyse the data from Nakagawa et al. (2015).
    -	[004_additional_sensitivity_analyses](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/blob/main/code/004_additional_sensitivity_analyses): used to run additional sensitivity analyses.
    -	[005_phylogeny_figures](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/blob/main/code/005_phylogeny_figures): used to generate Figure 1 (phylogenetic tree with effect sizes).

*	[data](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/tree/main/data): contains all [final](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/tree/main/data/final) datasets used in the analyses as well as the [phylogeny](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/tree/main/data/phylogeny) data generated and used for the phylogenetic multilevel meta-analyses. The folder [Nakagawa_et_al_2015](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/tree/main/data/Nakagawa_et_al_2015) contains the data corresponding to the re-analysis of the data provided by Nakagawa et al. (2015).

*	[figures](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/tree/main/figures): contains all figures exported (including the supplementary figures).

*	[models](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/tree/main/models): contains all models run to faciliate reproducibility.

## Metadata:
We provide the metadata for the raw dataset [here](https://github.com/ASanchez-Tojar/meta-analysis_canalisation_bird_wing/blob/main/data/METADATA.csv).

## Software and {Packages} used
R version 4.3.1 (2023-06-16 ucrt)

Platform: x86_64-w64-mingw32/x64 (64-bit)

Running under: Windows 10 x64 (build 19044)

Packages used: 
rsvg_2.6.1          
rgbif_3.8.1         
ggcorrplot_0.1.4.1  
clubSandwich_0.5.10 
orchaRd_2.0         
ggpubr_0.6.0       
RColorBrewer_1.1-3  
ggtree_3.10.1       
ggstance_0.3.6      
patchwork_1.3.0     
viridis_0.6.5       
viridisLite_0.4.2  
cowplot_1.1.1       
visdat_0.6.0        
ggsignif_0.6.4      
wesanderson_0.3.7   
diagram_1.6.5       
shape_1.4.6.1      
treebase_0.1.4      
metafor_4.6-0       
numDeriv_2016.8-1.1 
metadat_1.2-0       
Matrix_1.6-1        
ggstatsplot_0.12.5 
lubridate_1.9.4     
forcats_1.0.0       
stringr_1.5.1       
dplyr_1.1.4         
purrr_1.0.2         
readr_2.1.5        
tidyr_1.3.1         
tibble_3.2.1        
ggplot2_3.5.1       
tidyverse_2.0.0     
rotl_3.1.0          
ape_5.7-1          
