# Author's information: 
# Alfredo Sanchez-Tojar (@ASanchez_Tojar)
# Profile: https://goo.gl/PmpPEB
# Department of Evolutionary Biology, Bielefeld University (GER) 
# Email: alfredo.tojar@gmail.com
# Script first created on the 19th of July 2024

################################################################################
# Description of script and Instructions
################################################################################

# This script is for the analyses presented in section:

# S6. Re-analysis of Nakagawa’s et al. 2015 Example 2 dataset

# For the study:

# Klaus Reinhold, Alfredo Sánchez-Tójar. 2025. Wing length canalisation and 
# behaviour across birds: a phylogenetic meta-analysis of variance

################################################################################
# Packages needed
################################################################################

# remotes::install_github("thomasp85/patchwork")
# remotes::install_github("tidyverse/ggplot2", ref = remotes::github_pull("5592"))
# devtools::install_github("daniel1noble/orchaRd", force = TRUE)

pacman::p_load(ape,metafor,rotl, #openxlsx,readxl,fulltext,
               treebase,diagram,dplyr,stringr,ggplot2,wesanderson,
               ggsignif,visdat,cowplot,
               patchwork, #nice resourece page https://patchwork.data-imaginist.com/articles/guides/layout.html
               ggstance,ggtree,RColorBrewer,
               tidyverse,ggpubr,orchaRd,clubSandwich) #pander,brms,

# cleaning
rm(list=ls())


################################################################################
# Functions: 
################################################################################

# function to calculate heterogeneity from 
# https://github.com/Yefeng0920/heterogeneity_ecoevo/tree/main/function
# function to calculate heterogeneity - squared version, including CVH2 and M2
h.calc2 <- function(mod){
  # I2
  # sigma2_v = typical sampling error variance
  sigma2_v <- sum(1 / mod$vi) * (mod$k - 1) /
    (sum(1 / mod$vi)^2 - sum((1 / mod$vi)^2))
  # s^2_t = total variance
  I2_total <- 100 * (sum(mod$sigma2) / (sum(mod$sigma2) + sigma2_v))
  I2_each <- 100 * (mod$sigma2 / (sum(mod$sigma2) + sigma2_v))
  #names(I2_each) <- paste0("I2_", model$s.names)
  #names(I2_total) <- "I2_Total"
  I2s_Shinichi <- c(I2_total, I2_each)
  
  # matrix version  
  W <- solve(mod$V)
  X <- model.matrix(mod)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2_total2 <- 100* (sum(mod$sigma2) / (sum(mod$sigma2) + (mod$k - mod$p) / sum(diag(P))))
  I2_each2 <- 100* (mod$sigma2 / (sum(mod$sigma2) + (mod$k - mod$p) / sum(diag(P))))
  #names(I2_each2) <- paste0("I2_", model$s.names)
  #names(I2_total2) <- "I2_Total2"
  I2s_Wolfgang <- c(I2_total2, I2_each2)
  
  
  # CVH2
  CV_total <- (sum(mod$sigma2) / (mod$beta[1])^2)
  CV_each <- (mod$sigma2 / (mod$beta[1])^2)
  
  #names(CVB_each) <- paste0("CVB_", mod$s.names)
  #names(CVB_total) <- "CVB_total"
  CVHs <- c(CV_total, CV_each)
  
  # M2
  M_total <- (sum(mod$sigma2) / (sum(mod$sigma2) + (mod$beta[1])^2))
  M_each <- (mod$sigma2) / (sum(mod$sigma2) + (mod$beta[1])^2)
  Ms <- c(M_total, M_each)
  
  hs <- data.frame(I2s_Shinichi,CVHs,Ms)
  rownames(hs) <- c("Total", mod$s.names)
  return(hs)
}

# function to estimate typical sampling error variance from 
# https://github.com/Yefeng0920/heterogeneity_ecoevo/tree/main/function
sigma2_v <- function(mod){
  sigma2_v <- sum(1 / mod$vi) * (mod$k - 1) /
    (sum(1 / mod$vi)^2 - sum((1 / mod$vi)^2))
  return(sigma2_v)
}


# function for extracting mean and CI from each metafor model
estimates.CI <- function(model){
  db.mf <- data.frame(model$b,row.names = 1:nrow(model$b))
  db.mf <- cbind(db.mf,model$ci.lb,model$ci.ub,row.names(model$b),model$pval)
  names(db.mf) <- c("mean","lower","upper","estimate","p-value")
  return(db.mf[,c("estimate","mean","lower","upper","p-value")])
}

################################################################################
# Importing the data 
################################################################################

# importing data
body.variance.full <- read.csv("data/Nakagawa_et_al_2015/mee312309-sup-0010-datas10.csv",
                               header=T,sep=",")

# reformatting some variables types to factors
cols.factors <- c("species","animal","reference","unit","Trait","trait_details")
body.variance.full[cols.factors] <- lapply(body.variance.full[cols.factors], 
                                           factor)


# species
table(sort(body.variance.full$species))
sort(table(body.variance.full$species))
length(unique(body.variance.full$species))
summary(as.numeric(table(body.variance.full$species)))
sd(as.numeric(table(body.variance.full$species)))

# trait
table(body.variance.full$trait)
table(body.variance.full$trait_details)

body.variance.full$trait.updated <-
  recode(body.variance.full$trait,
         "Fl\xf4gell_nge" = "wing length",
         "body length" = "body size",
         "body weight" = "body mass",
         "Body weight" = "body mass",
         "culemn length" = "beak size",
         "Culmen" = "beak size",
         "culmen length" = "beak size",
         "Culmen length" = "beak size",
         "flipper length" = "wing length",
         "head length" = "head size",
         "head width" = "head size",
         "mass" = "body mass",
         "skull length" = "head size",
         "total culmen " = "beak size",
         "total length" = "body size",
         "weight" = "body mass",
         "weigth " = "body mass",
         "wing" = "wing length",
         "Wing" = "wing length",
         "wing " = "wing length",
         "wing chord" = "wing length",
         "Wing length" = "wing length",
         .default = levels(body.variance.full$trait))

table(body.variance.full$trait.updated)

head(body.variance.full)


body.variance.full$pairID <- 1:nrow(body.variance.full)

# splitting the database so that males and females can be analyzed together
body.variance.full.final.males <- as.data.frame(body.variance.full 
                                                %>% select(-female_mean, 
                                                           -female_n, 
                                                           -female_SD,
                                                           -X) 
                                                %>% filter(!is.na(male_mean))#excluding NA values (including those assigned above)
                                                %>% rename(mean = male_mean,
                                                           n = male_n, 
                                                           SD = male_SD)
                                                %>% mutate(sex = "male"))

nrow(body.variance.full.final.males)

body.variance.full.final.females <- as.data.frame(body.variance.full 
                                                  %>% select(-male_mean, 
                                                             -male_n, 
                                                             -male_SD,
                                                             -X)
                                                  %>% filter(!is.na(female_mean))#excluding NA values (including those assigned above)
                                                  %>% rename(mean = female_mean,
                                                             n = female_n,
                                                             SD = female_SD)
                                                  %>% mutate(sex = "female"))


nrow(body.variance.full.final.females)



# building the combined male-female database for the analyses
body.variance.full.sex.sep <- rbind(body.variance.full.final.males,
                                    body.variance.full.final.females)
body.variance.full.sex.sep$sex <- as.factor(body.variance.full.sex.sep$sex)

summary(body.variance.full.sex.sep)

# How many studies in total?
length(unique(body.variance.full.sex.sep$reference))


# obtaining dataframe listing the Open Tree identifiers potentially matching
# our list of species (be aware that this will take a few minutes. You can avoid
# waiting by loading the data below)

body.variance.full.sex.sep$species <- str_replace(body.variance.full.sex.sep$species,"_"," ")
study.species <- as.character(unique(body.variance.full.sex.sep$species))

# How many species?
length(unique(study.species))

# taxa <- tnrs_match_names(names = study.species)
# 
# # check approximate matches
# taxa[taxa$approximate_match==TRUE,]
# 
# # check approximate matches
# taxa[taxa$is_synonym==TRUE,]
# 
# # check number of matches
# taxa[taxa$number_matches>1,]

# Based on the checks above, we had to implement some changes to the species names:
# The remaining changes are typos that needed correction, update in species names or bringing subspecies up to species level so that everything is comparable
body.variance.full.sex.sep$species.updated <-
  recode(body.variance.full.sex.sep$species,
         "Puffinus iherninieri" = "Puffinus lherminieri",
         "Phalacocorax atriceps" = "Leucocarbo atriceps",
         "Calidris marittima" = "Calidris maritima",
         "Calidris temmickii" = "Calidris temminckii",
         "Accipter cooperii" = "Accipiter cooperii",
         "Dendrocopus medius" = "Dendrocopos medius",
         "Phalacrocorax aristotelis" = "Gulosus aristotelis",
         "Pygoscelis antarctica" = "Pygoscelis antarcticus",
         "Dendrocopos pubescens" = "Picoides pubescens",
         "Phylloscopus inornatus" = "Abrornis inornata",
         "Parus caeruleus" = "Cyanistes caeruleus",
         "Carduelis chloris" = "Chloris chloris",
         "Cinclorhamphus cruralis" = "Cincloramphus cruralis",
         "Puffinus carneipes" = "Ardenna carneipes",
         "Puffinus creatopus" = "Ardenna creatopus",
         "Puffinus gravis" = "Ardenna gravis",
         "Puffinus tenuirostris" = "Ardenna tenuirostris",
         "Puffinus bulleri" = "Ardenna bulleri",
         "Puffinus griseus" = "Ardenna grisea",
         "Puffinus pacificus" = "Ardenna pacifica",
         "Bonasa sewerzowi" = "Tetrastes sewerzowi",
         "Carpodacus mexicanus" = "Haemorhous mexicanus",
         "Catharacta skua" = "Stercorarius skua",
         "Diomedea nigripes" = "Phoebastria nigripes",
         "Diomedea immutabilis" = "Phoebastria immutabilis",
         "Miliaria calandra" = "Emberiza calandra",
         .default = levels(body.variance.full.sex.sep$species))

body.variance.full.sex.sep$species.updated <- factor(body.variance.full.sex.sep$species.updated)
length(unique(body.variance.full.sex.sep$species.updated))



# repeating the process after fixing the species names
# obtaining dataframe listing the Open Tree identifiers potentially matching our
# list of species (be aware that this will take a few minutes. You can avoid
# waiting by loading the data below)
study.species <- as.character(unique(body.variance.full.sex.sep$species.updated))

# taxa <- tnrs_match_names(names = study.species)
# 
# # check approximate matches
# taxa[taxa$approximate_match==TRUE,]
# 
# # check approximate matches
# taxa[taxa$is_synonym==TRUE,]
# 
# # check number of matches
# taxa[taxa$number_matches>1,]
# 
# # some further checks
# ott_id_tocheck <- taxa[taxa$number_matches != 1,"ott_id"]
# 
# for(i in 1:length(ott_id_tocheck)){
#   print(inspect(taxa, ott_id = ott_id_tocheck[i]))
# }
# 
# # saving the taxonomic data last created on ... to speed the process in the
# # future and allow full computional reproducibility
# save(taxa,file = "data/Nakagawa_et_al_2015/taxa_Open_Tree_of_Life_Nakagawa_et_al_2015_20240719.RData")

# loading the most updated taxonomic data created
load("data/Nakagawa_et_al_2015/taxa_Open_Tree_of_Life_Nakagawa_et_al_2015_20240719.RData") #taxa


# everything seems in order now

# # retrieving phylogenetic relationships among taxa in the form of a trimmed
# # sub-tree
# tree <- tol_induced_subtree(ott_ids = taxa[["ott_id"]], label_format = "name")
# 
# # we need to check for the existence of polytomies
# is.binary(tree)
# 
# # to confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# 
# tree$tip.label <- gsub("_"," ", tree$tip.label)
# intersect(as.character(tree$tip.label),
#           as.character(body.variance.full.sex.sep$species.updated))
# setdiff(as.character(body.variance.full.sex.sep$species.updated),
#         as.character(tree$tip.label)) #listed in our database but not in the tree
# setdiff(as.character(tree$tip.label),
#         as.character(body.variance.full.sex.sep$species.updated)) # listed in the tree but not in our database
# 
# 
# # we can now save the tree
# save(tree, file = "data/Nakagawa_et_al_2015/tree_random_Nakagawa_et_al_2015_20240719.RData")
# 
# 
# # compute branch lengths of tree
# phylo_branch <- compute.brlen(tree, method = "Grafen", power = 1)
# 
# # check tree is ultrametric
# is.ultrametric(phylo_branch) # TRUE
# 
# # matrix to be included in the models
# phylo_cor <- vcv(phylo_branch, cor = T)
# 
# # finally, save matrix for future analyses
# save(phylo_cor, file = "data/Nakagawa_et_al_2015/phylo_cor_Nakagawa_et_al_2015_20240719.Rdata")


# loading tree and matrix
load("data/Nakagawa_et_al_2015/tree_random_Nakagawa_et_al_2015_20240719.RData")
load("data/Nakagawa_et_al_2015/phylo_cor_Nakagawa_et_al_2015_20240719.Rdata")



################################################################################
# Calculating lnCV
################################################################################

# calculating lnCV
body.variance.full.sex.sep <- as.data.frame(escalc(measure = "CVLN", 
                                                   mi=mean, 
                                                   sdi=SD, 
                                                   ni=n,
                                                   data=body.variance.full.sex.sep,
                                                   add.measure=FALSE,
                                                   append=TRUE))

summary(body.variance.full.sex.sep)

################################################################################
# Meta-analysis
################################################################################

# setting reference levels for the factorial moderators
body.variance.full.sex.sep$trait.updated <- relevel(body.variance.full.sex.sep$trait.updated, 
                                                    ref = "wing length")

table(body.variance.full.sex.sep$trait.updated)


# generating a unit level effect to model residual/within-study variance
body.variance.full.sex.sep$effectsizeID <- 1:nrow(body.variance.full.sex.sep)

# generating a new species variable to account for repeated values for some 
# species + phylogeny (sensu Cinar et al. 2020, MEE)
body.variance.full.sex.sep$species.updated.rep <- body.variance.full.sex.sep$species.updated


VCV.0.50.study_ID <- impute_covariance_matrix(vi = body.variance.full.sex.sep$vi, 
                                              cluster = body.variance.full.sex.sep$reference, 
                                              r = 0.50)


################################################################################
# intercept-only but with VCV 0.5 for study_ID
# model_intercept_only_Nakawaga_VCV050_studyID <- rma.mv(yi=yi,
#                                                        V = VCV.0.50.study_ID,
#                                                        random = list(~ 1 | species.updated.rep,
#                                                                      ~ 1 | species.updated,
#                                                                      ~ 1 | reference,
#                                                                      ~ 1 | pairID,
#                                                                      ~ 1 | effectsizeID),
#                                                        R = list(species.updated = phylo_cor),
#                                                        data=body.variance.full.sex.sep,
#                                                        method="REML",
#                                                        test="t")
# 
# save(model_intercept_only_Nakawaga_VCV050_studyID,
#      file="models/model_intercept_only_Nakawaga_VCV050_studyID.RData")
load("models/model_intercept_only_Nakawaga_VCV050_studyID.RData")

summary(model_intercept_only_Nakawaga_VCV050_studyID,3)
predict(model_intercept_only_Nakawaga_VCV050_studyID, transf=exp)


################################################################################
# HETEROGENEITY for our model of choice
################################################################################

# # I2, CV and M (Yang et al. 2023)
# round(h.calc(model_intercept_only_Nakawaga_VCV050_studyID),3)

# same but using orchaRd, where m1_ml()'s bug has now been fixed
round(i2_ml(model_intercept_only_Nakawaga_VCV050_studyID),1)
round(cvh2_ml(model_intercept_only_Nakawaga_VCV050_studyID),3)
round(m2_ml(model_intercept_only_Nakawaga_VCV050_studyID),3)

# total variance
round(sum(model_intercept_only_Nakawaga_VCV050_studyID$sigma2),3)

# Typical sampling variance (which captures the statistical noise of the data, 
# but is rarely reported in the current meta-analytic practice):
sigma2_v(model_intercept_only_Nakawaga_VCV050_studyID)


################################################################################
# Visualize heterogeneity
## make dataframe
h_status <- h.calc2(model_intercept_only_Nakawaga_VCV050_studyID)

# adding sigmas
h_status$sigma2s <- c(sum(model_intercept_only_Nakawaga_VCV050_studyID$sigma2),
                      model_intercept_only_Nakawaga_VCV050_studyID$sigma2[1],
                      model_intercept_only_Nakawaga_VCV050_studyID$sigma2[2],
                      model_intercept_only_Nakawaga_VCV050_studyID$sigma2[3],
                      model_intercept_only_Nakawaga_VCV050_studyID$sigma2[4],
                      model_intercept_only_Nakawaga_VCV050_studyID$sigma2[5])

round(h_status$I2s_Shinichi,1)
round(h_status$CVHs,3)
round(h_status$Ms,3)
round(h_status$sigma2s,3)

h_status$levels <- rownames(h_status)  
h_status$levels <- dplyr::recode(h_status$levels, 
                                 "Total" = "Total",  
                                 "species.updated.rep" = "Species",  
                                 "species.updated" = "Phylogeny", 
                                 "reference" = "Between-study", 
                                 "pairID" = "Paired estimates",
                                 "effectsizeID" = "Within-study")
h_status$levels <- as.factor(h_status$levels)
h_status$levels <- factor(h_status$levels, levels = c("Total",
                                                      "Between-study",
                                                      "Paired estimates",
                                                      "Phylogeny", 
                                                      "Species", 
                                                      "Within-study"))


p.sigma <- ggplot(h_status, aes(levels, sigma2s)) +
  geom_col(alpha = 1,
           color = wes_palette('GrandBudapest1', 4, type = 'discrete')[1],
           fill = wes_palette('GrandBudapest1', 4, type = 'discrete')[1]) +
  labs(y = expression("Variance"), x = "" ,
       title = "Unstandardised heterogeneity metric") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + # https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
  theme_bw() +
  theme(legend.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")
  )

# I2
p.I2 <- ggplot(h_status, aes(levels, I2s_Shinichi/100)) +
  geom_col(alpha = 1, color = wes_palette('GrandBudapest1', 4, 
                                          type = 'discrete')[2], 
           fill = wes_palette('GrandBudapest1', 4, type = 'discrete')[2]) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = expression(paste(italic(I)[]^2)), x = "" , 
       title = "Source of heterogeneity") + 
  theme_bw() +
  theme(legend.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")
  ) 

# CV
p.CV <- ggplot(h_status, aes(levels, CVHs)) +
  geom_col(alpha = 1, color = wes_palette('GrandBudapest1', 4, 
                                          type = 'discrete')[3], 
           fill = wes_palette('GrandBudapest1', 4, type = 'discrete')[3]) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  
  labs(y = expression(paste(italic(CVH2)[])), x = "Strata" , 
       title = "Magnitude of heterogeneity") + 
  theme_bw() +
  theme(legend.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")
  )   


# M
p.M <- ggplot(h_status, aes(x = levels, y = Ms)) +
  geom_col(alpha = 1, color = wes_palette('GrandBudapest1', 4, 
                                          type = 'discrete')[4], 
           fill = wes_palette('GrandBudapest1', 4, type = 'discrete')[4]) +
  labs(y = expression(paste(italic(M2)[])), x = "Strata" , 
       title = "Magnitude of heterogeneity") + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + # https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
  theme_bw() +
  theme(legend.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")
  )  


# saving the figure
png(filename = 'figures/supplementary_figures/Figure_S12_Heterogeneity_pluralistic_approach.png', 
    width = 19, height = 19, units = 'cm', 
    res = 600)

#png(filename = "all.png", width = 9, height = 8, units = "in", type = "windows", res = 400)
p.sigma + p.I2 + p.CV + p.M + plot_layout(ncol = 2, nrow = 2) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold"),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1))

dev.off()


################################################################################
# Meta-regressions
################################################################################

################################################################################
# Confirming sex effect (females larger CV than males) from: 
# Reinhold & Engqvist 2013: https://doi.org/10.1111/evo.12224 (although multiple traits included not only wing size)
# Nakagawa et al. 2015: https://doi.org/10.1111/2041-210X.12309 (although multiple traits included not only wing size)

# metaregression_unimoderator_sex_Nakagawa <- rma.mv(yi = yi,
#                                                    V = VCV.0.50.study_ID,
#                                                    mod = ~ 1 + sex,
#                                                    random = list(~ 1 | species.updated.rep,
#                                                                  ~ 1 | species.updated,
#                                                                  ~ 1 | reference,
#                                                                  ~ 1 | pairID,
#                                                                  ~ 1 | effectsizeID),
#                                                    R = list(species.updated = phylo_cor),
#                                                    data=body.variance.full.sex.sep,
#                                                    method="REML",
#                                                    test="t")
# 
# save(metaregression_unimoderator_sex_Nakagawa,
#      file="models/metaregression_unimoderator_sex_Nakagawa.RData")
load("models/metaregression_unimoderator_sex_Nakagawa.RData")

summary(metaregression_unimoderator_sex_Nakagawa)

# getting marginal R2
round(r2_ml(metaregression_unimoderator_sex_Nakagawa)*100,2)


# the equivalent model but without intercept
# metaregression_unimoderator_sex_Nakagawa.nointercept <- rma.mv(yi = yi,
#                                                                V = VCV.0.50.study_ID,
#                                                                mod = ~ -1 + sex,
#                                                                random = list(~ 1 | species.updated.rep,
#                                                                              ~ 1 | species.updated,
#                                                                              ~ 1 | reference,
#                                                                              ~ 1 | pairID,
#                                                                              ~ 1 | effectsizeID),
#                                                                R = list(species.updated = phylo_cor),
#                                                                data=body.variance.full.sex.sep,
#                                                                method="REML",
#                                                                test="t")
# 
# save(metaregression_unimoderator_sex_Nakagawa.nointercept,
#      file="models/metaregression_unimoderator_sex_Nakagawa_nointercept.RData")
load("models/metaregression_unimoderator_sex_Nakagawa_nointercept.RData")

summary(metaregression_unimoderator_sex_Nakagawa.nointercept)
round(exp(metaregression_unimoderator_sex_Nakagawa.nointercept$beta)*100,2)
round(exp(metaregression_unimoderator_sex_Nakagawa.nointercept$ci.lb)*100,2)
round(exp(metaregression_unimoderator_sex_Nakagawa.nointercept$ci.ub)*100,2)


summary(body.variance.full.sex.sep$yi)
exp(min(body.variance.full.sex.sep$yi))*100
exp(max(body.variance.full.sex.sep$yi))*100

summary(attributes(body.variance.full.sex.sep$yi)$ni)



################################################################################
# plotting results

orchard.plot.sex <- orchaRd::orchard_plot(metaregression_unimoderator_sex_Nakagawa, 
                                          mod = "sex", 
                                          #group = "species.updated.rep", 
                                          group = "trait.updated", 
                                          xlab = "Effect size (CV %)",
                                          trunk.size = 2,
                                          branch.size = 3.5,
                                          twig.size = 1.25,
                                          g = FALSE,
                                          transfm  = "percent",
                                          colour = TRUE)

# saving the figure
png(filename = 'figures/supplementary_figures/Figure_S13_Sex_Chromosome_Hypothesis_Nakagawa_et_al_2015_overall.png',
    width = 21, height = 10, units = 'cm',
    res = 600)

orchard.plot.sex + 
  scale_y_continuous(limits = c(0,31),
                     breaks = seq(0,31,1)) + theme()

dev.off()




################################################################################
# repeating all analyses with wing length only
################################################################################

body.variance.full.sex.sep.wing.length <- body.variance.full.sex.sep[body.variance.full.sex.sep$trait.updated=="wing length",]



# repeating the process after fixing the species names
# obtaining dataframe listing the Open Tree identifiers potentially matching our
# list of species (be aware that this will take a few minutes. You can avoid
# waiting by loading the data below)
study.species <- as.character(unique(body.variance.full.sex.sep.wing.length$species.updated))

# taxa <- tnrs_match_names(names = study.species)
# 
# # check approximate matches
# taxa[taxa$approximate_match==TRUE,]
# 
# # check approximate matches
# taxa[taxa$is_synonym==TRUE,]
# 
# # check number of matches
# taxa[taxa$number_matches>1,]
# 
# # some further checks
# ott_id_tocheck <- taxa[taxa$number_matches != 1,"ott_id"]
# 
# for(i in 1:length(ott_id_tocheck)){
#   print(inspect(taxa, ott_id = ott_id_tocheck[i]))
# }
# 
# # saving the taxonomic data last created on ... to speed the process in the
# # future and allow full computional reproducibility
# save(taxa,file = "data/Nakagawa_et_al_2015/taxa_Open_Tree_of_Life_Nakagawa_et_al_2015_wing_length_20240719.RData")

# loading the most updated taxonomic data created
load("data/Nakagawa_et_al_2015/taxa_Open_Tree_of_Life_Nakagawa_et_al_2015_wing_length_20240719.RData") #taxa


# everything seems in order now

# # retrieving phylogenetic relationships among taxa in the form of a trimmed
# # sub-tree
# tree <- tol_induced_subtree(ott_ids = taxa[["ott_id"]], label_format = "name")
# 
# # we need to check for the existence of polytomies
# is.binary(tree)
# 
# # to confirm that our tree covers all the species we wanted it to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# 
# tree$tip.label <- gsub("_"," ", tree$tip.label)
# intersect(as.character(tree$tip.label),
#           as.character(body.variance.full.sex.sep.wing.length$species.updated))
# setdiff(as.character(body.variance.full.sex.sep.wing.length$species.updated),
#         as.character(tree$tip.label)) #listed in our database but not in the tree
# setdiff(as.character(tree$tip.label),
#         as.character(body.variance.full.sex.sep.wing.length$species.updated)) # listed in the tree but not in our database
# 
# 
# # we can now save the tree
# save(tree, file = "data/Nakagawa_et_al_2015/tree_random_Nakagawa_et_al_2015_wing_length_20240719.RData")
# 
# 
# # compute branch lengths of tree
# phylo_branch <- compute.brlen(tree, method = "Grafen", power = 1)
# 
# # check tree is ultrametric
# is.ultrametric(phylo_branch) # TRUE
# 
# # matrix to be included in the models
# phylo_cor <- vcv(phylo_branch, cor = T)
# 
# # finally, save matrix for future analyses
# save(phylo_cor, file = "data/Nakagawa_et_al_2015/phylo_cor_Nakagawa_et_al_2015_wing_length_20240719.Rdata")


# loading tree and matrix
load("data/Nakagawa_et_al_2015/tree_random_Nakagawa_et_al_2015_wing_length_20240719.RData")
load("data/Nakagawa_et_al_2015/phylo_cor_Nakagawa_et_al_2015_wing_length_20240719.Rdata")



################################################################################
# Calculating lnCV
################################################################################

# calculating lnCV
body.variance.full.sex.sep.wing.length <- as.data.frame(escalc(measure = "CVLN", 
                                                               mi=mean, 
                                                               sdi=SD, 
                                                               ni=n,
                                                               data=body.variance.full.sex.sep.wing.length,
                                                               add.measure=FALSE,
                                                               append=TRUE))

summary(body.variance.full.sex.sep.wing.length)
length(unique(body.variance.full.sex.sep.wing.length$reference))
length(unique(body.variance.full.sex.sep.wing.length$species.updated))

################################################################################
# Meta-analysis
################################################################################

VCV.0.50.study_ID <- impute_covariance_matrix(vi = body.variance.full.sex.sep.wing.length$vi, 
                                              cluster = body.variance.full.sex.sep.wing.length$reference, 
                                              r = 0.50)


################################################################################
# intercept-only but with VCV 0.5 for study_ID
# model_intercept_only_Nakawaga_wing_length_VCV050_studyID <- rma.mv(yi=yi,
#                                                                    V = VCV.0.50.study_ID,
#                                                                    random = list(~ 1 | species.updated.rep,
#                                                                                  ~ 1 | species.updated,
#                                                                                  ~ 1 | reference,
#                                                                                  ~ 1 | pairID,
#                                                                                  ~ 1 | effectsizeID),
#                                                                    R = list(species.updated = phylo_cor),
#                                                                    data=body.variance.full.sex.sep.wing.length,
#                                                                    method="REML",
#                                                                    test="t")
# 
# save(model_intercept_only_Nakawaga_wing_length_VCV050_studyID,
#      file="models/model_intercept_only_Nakawaga_wing_length_VCV050_studyID.RData")
load("models/model_intercept_only_Nakawaga_wing_length_VCV050_studyID.RData")

summary(model_intercept_only_Nakawaga_wing_length_VCV050_studyID,3)
predict(model_intercept_only_Nakawaga_wing_length_VCV050_studyID, transf=exp)


################################################################################
# HETEROGENEITY for our model of choice
################################################################################

# # I2, CV and M (Yang et al. 2023)
# round(h.calc(model_intercept_only_Nakawaga_VCV050_studyID),3)

# same but using orchaRd, where m1_ml()'s bug has now been fixed
round(i2_ml(model_intercept_only_Nakawaga_wing_length_VCV050_studyID),1)
round(cvh2_ml(model_intercept_only_Nakawaga_wing_length_VCV050_studyID),3)
round(m2_ml(model_intercept_only_Nakawaga_wing_length_VCV050_studyID),3)

# total variance
round(sum(model_intercept_only_Nakawaga_wing_length_VCV050_studyID$sigma2),3)

# Typical sampling variance (which captures the statistical noise of the data, 
# but is rarely reported in the current meta-analytic practice):
sigma2_v(model_intercept_only_Nakawaga_wing_length_VCV050_studyID)


################################################################################
# Visualize heterogeneity
## make dataframe
h_status <- h.calc2(model_intercept_only_Nakawaga_wing_length_VCV050_studyID)

# adding sigmas
h_status$sigma2s <- c(sum(model_intercept_only_Nakawaga_wing_length_VCV050_studyID$sigma2),
                      model_intercept_only_Nakawaga_wing_length_VCV050_studyID$sigma2[1],
                      model_intercept_only_Nakawaga_wing_length_VCV050_studyID$sigma2[2],
                      model_intercept_only_Nakawaga_wing_length_VCV050_studyID$sigma2[3],
                      model_intercept_only_Nakawaga_wing_length_VCV050_studyID$sigma2[4],
                      model_intercept_only_Nakawaga_wing_length_VCV050_studyID$sigma2[5])

round(h_status$I2s_Shinichi,1)
round(h_status$CVHs,3)
round(h_status$Ms,3)
round(h_status$sigma2s,3)

h_status$levels <- rownames(h_status)  
h_status$levels <- dplyr::recode(h_status$levels, 
                                 "Total" = "Total",  
                                 "species.updated.rep" = "Species",  
                                 "species.updated" = "Phylogeny", 
                                 "reference" = "Between-study", 
                                 "pairID" = "Paired estimates",
                                 "effectsizeID" = "Within-study")
h_status$levels <- as.factor(h_status$levels)
h_status$levels <- factor(h_status$levels, levels = c("Total",
                                                      "Between-study",
                                                      "Paired estimates",
                                                      "Phylogeny", 
                                                      "Species", 
                                                      "Within-study"))


p.sigma <- ggplot(h_status, aes(levels, sigma2s)) +
  geom_col(alpha = 1,
           color = wes_palette('GrandBudapest1', 4, type = 'discrete')[1],
           fill = wes_palette('GrandBudapest1', 4, type = 'discrete')[1]) +
  labs(y = expression("Variance"), x = "" ,
       title = "Unstandardised heterogeneity metric") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + # https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
  theme_bw() +
  theme(legend.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")
  )

# I2
p.I2 <- ggplot(h_status, aes(levels, I2s_Shinichi/100)) +
  geom_col(alpha = 1, color = wes_palette('GrandBudapest1', 4, 
                                          type = 'discrete')[2], 
           fill = wes_palette('GrandBudapest1', 4, type = 'discrete')[2]) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = expression(paste(italic(I)[]^2)), x = "" , 
       title = "Source of heterogeneity") + 
  theme_bw() +
  theme(legend.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")
  ) 

# CV
p.CV <- ggplot(h_status, aes(levels, CVHs)) +
  geom_col(alpha = 1, color = wes_palette('GrandBudapest1', 4, 
                                          type = 'discrete')[3], 
           fill = wes_palette('GrandBudapest1', 4, type = 'discrete')[3]) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) + # https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
  
  labs(y = expression(paste(italic(CVH2)[])), x = "Strata" , 
       title = "Magnitude of heterogeneity") + 
  theme_bw() +
  theme(legend.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")
  )   


# M
p.M <- ggplot(h_status, aes(x = levels, y = Ms)) +
  geom_col(alpha = 1, color = wes_palette('GrandBudapest1', 4, 
                                          type = 'discrete')[4], 
           fill = wes_palette('GrandBudapest1', 4, type = 'discrete')[4]) +
  labs(y = expression(paste(italic(M2)[])), x = "Strata" , 
       title = "Magnitude of heterogeneity") + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) + # https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
  theme_bw() +
  theme(legend.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black")
  )  


# saving the figure
png(filename = 'figures/supplementary_figures/Figure_S14_Heterogeneity_pluralistic_approach.png', 
    width = 19, height = 19, units = 'cm', 
    res = 600)

#png(filename = "all.png", width = 9, height = 8, units = "in", type = "windows", res = 400)
p.sigma + p.I2 + p.CV + p.M + plot_layout(ncol = 2, nrow = 2) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold"),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1))

dev.off()



################################################################################
# Meta-regressions
################################################################################

################################################################################
# Confirming sex effect (females larger CV than males) from: 
# Reinhold & Engqvist 2013: https://doi.org/10.1111/evo.12224 (although multiple traits included not only wing size)
# Nakagawa et al. 2015: https://doi.org/10.1111/2041-210X.12309 (although multiple traits included not only wing size)

# metaregression_unimoderator_sex_Nakagawa_wing_length <- rma.mv(yi = yi,
#                                                    V = VCV.0.50.study_ID,
#                                                    mod = ~ 1 + sex,
#                                                    random = list(~ 1 | species.updated.rep,
#                                                                  ~ 1 | species.updated,
#                                                                  ~ 1 | reference,
#                                                                  ~ 1 | pairID,
#                                                                  ~ 1 | effectsizeID),
#                                                    R = list(species.updated = phylo_cor),
#                                                    data=body.variance.full.sex.sep.wing.length,
#                                                    method="REML",
#                                                    test="t")
# 
# save(metaregression_unimoderator_sex_Nakagawa_wing_length,
#      file="models/metaregression_unimoderator_sex_Nakagawa_wing_length.RData")
load("models/metaregression_unimoderator_sex_Nakagawa_wing_length.RData")

summary(metaregression_unimoderator_sex_Nakagawa_wing_length)

# getting marginal R2
round(r2_ml(metaregression_unimoderator_sex_Nakagawa_wing_length)*100,2)


# # the equivalent model but without intercept
# metaregression_unimoderator_sex_Nakagawa.nointercept_wing_length <- rma.mv(yi = yi,
#                                                                V = VCV.0.50.study_ID,
#                                                                mod = ~ -1 + sex,
#                                                                random = list(~ 1 | species.updated.rep,
#                                                                              ~ 1 | species.updated,
#                                                                              ~ 1 | reference,
#                                                                              ~ 1 | pairID,
#                                                                              ~ 1 | effectsizeID),
#                                                                R = list(species.updated = phylo_cor),
#                                                                data=body.variance.full.sex.sep.wing.length,
#                                                                method="REML",
#                                                                test="t")
# 
# save(metaregression_unimoderator_sex_Nakagawa.nointercept_wing_length,
#      file="models/metaregression_unimoderator_sex_Nakagawa_nointercept_wing_length.RData")
load("models/metaregression_unimoderator_sex_Nakagawa_nointercept_wing_length.RData")

summary(metaregression_unimoderator_sex_Nakagawa.nointercept_wing_length)
round(exp(metaregression_unimoderator_sex_Nakagawa.nointercept_wing_length$beta)*100,2)
round(exp(metaregression_unimoderator_sex_Nakagawa.nointercept_wing_length$ci.lb)*100,2)
round(exp(metaregression_unimoderator_sex_Nakagawa.nointercept_wing_length$ci.ub)*100,2)


summary(body.variance.full.sex.sep.wing.length$yi)
exp(min(body.variance.full.sex.sep.wing.length$yi))*100
exp(max(body.variance.full.sex.sep.wing.length$yi))*100

summary(attributes(body.variance.full.sex.sep.wing.length$yi)$ni)



################################################################################
# plotting results

orchard.plot.sex <- orchaRd::orchard_plot(metaregression_unimoderator_sex_Nakagawa_wing_length, 
                                          mod = "sex", 
                                          group = "species.updated.rep", 
                                          xlab = "Effect size (CV %)",
                                          trunk.size = 2,
                                          branch.size = 3.5,
                                          twig.size = 1.25,
                                          transfm  = "percent",
                                          colour = F)

# saving the figure
png(filename = 'figures/supplementary_figures/Figure_S15_Sex_Chromosome_Hypothesis_Nakagawa_et_al_2015_wing_length.png',
    width = 21, height = 10, units = 'cm',
    res = 600)

orchard.plot.sex + 
  scale_y_continuous(limits = c(0,10),
                     breaks = seq(0,10,1))

dev.off()
