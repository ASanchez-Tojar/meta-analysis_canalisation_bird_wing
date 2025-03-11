################################################################################
# Author's information: 
# Alfredo Sanchez-Tojar (@ASanchez_Tojar)
# Profile: https://goo.gl/PmpPEB
# Department of Evolutionary Biology, Bielefeld University (GER) 
# Email: alfredo.tojar@gmail.com
# Script first created on the 25th of August 2021

################################################################################
# Description of script and Instructions
################################################################################

# This script imports data on bird wing length to clean it and get it read it
# for calculating and analysing wing length coefficient of variations across 
# bird species, and understand if differences across species can be explained by 
# species-specific and population-specific differences in migratory status and 
# feeding type.

# For the study:

# Klaus Reinhold, Alfredo Sánchez-Tójar. 2025. Wing length canalisation and 
# behaviour across birds: a phylogenetic meta-analysis of variance

################################################################################
# Packages needed
################################################################################

# install.packages("pacman")
pacman::p_load(ape,rotl,tidyverse,ggstatsplot)

# cleaning environment
rm(list=ls())


################################################################################
# Functions: 
################################################################################

# none needed

################################################################################
# Importing the data 
################################################################################

# importing final dataset containing means, SDs, Ns from each species and pop
body.variance.full <- read.csv("data/final/01_final_bird_wing_size_dataset.csv",
                               header=T,sep=",")

# converting some variables types to factors
cols.factors <- c("species","source","reference_link","population_ID",
                  "reference","trait","unit","exclude","location","Measurement")
body.variance.full[cols.factors] <- lapply(body.variance.full[cols.factors], 
                                           factor)

################################################################################
# Cleaning up the dataset: variable by variable
################################################################################

################################################################################
# species
table(sort(body.variance.full$species))
sort(table(body.variance.full$species))
length(unique(body.variance.full$species)) # note that subspecies are counted
summary(as.numeric(table(body.variance.full$species)))
sd(as.numeric(table(body.variance.full$species)))


################################################################################
# trait
table(body.variance.full$trait)
# in reality they all mean more or less the same. Flattened vs unflattened is 
# the main difference and we won't make a difference for the analyses.


################################################################################
# mean.male and mean.female
summary(body.variance.full$mean.male)
summary(body.variance.full$mean.female)

# # checking for gross outliers
# plot(body.variance.full$mean.male,body.variance.full$mean.female,
#      main = "Wing length Female~Male Overall correlation",
#      xlab = "Mean Wing length (Males)",
#      ylab = "Mean Wing length (Females)",
#      xlim = c(0,650),
#      ylim = c(0,650))
# abline(1,1)

# only bustards are a bit of an outlier, but it makes sense
# body.variance.full[body.variance.full$mean.male>550,]

# calculating the correlation value
cor.test(body.variance.full$mean.male,body.variance.full$mean.female)


################################################################################
# n.male
summary(body.variance.full$n.male)
table(!(is.na(body.variance.full$n.male)))

# n.female
summary(body.variance.full$n.female)
table(!(is.na(body.variance.full$n.female)))

# "extreme" values for Euplectes orix, from Craig et al. 2010; and also for 
# Ficedula albicollis: double-checking the original references indicates that 
# they are correct
body.variance.full[body.variance.full$n.male>3000 & !(is.na(body.variance.full$n.male)),]
body.variance.full[body.variance.full$n.female>3000&  !(is.na(body.variance.full$n.female)),]


################################################################################
# SD.male
summary(body.variance.full$SD.male)
table(!(is.na(body.variance.full$SD.male)))

# SD.female
summary(body.variance.full$SD.female)
table(!(is.na(body.variance.full$SD.female)))

# checking an SD == 0 that is suspicious
body.variance.full[body.variance.full$SD.male==0 & !(is.na(body.variance.full$SD.male)),] 
# sample size is very small, so by chance, the three individuals could have had 
# indeed the same exact measurement. We will exclude it from our analyses.

body.variance.full[body.variance.full$SD.female==0 & !(is.na(body.variance.full$SD.female)),] 


# assigning NA to mean, SD and n to any estimate with an n == 1
body.variance.full$mean.male <- ifelse(body.variance.full$n.male==1,
                                       NA,
                                       body.variance.full$mean.male)

body.variance.full$n.male <- ifelse(body.variance.full$n.male==1,
                                    NA,
                                    body.variance.full$n.male)

body.variance.full$SD.male <- ifelse(body.variance.full$n.male==1,
                                     NA,
                                     body.variance.full$SD.male)

body.variance.full$mean.female <- ifelse(body.variance.full$n.female==1,
                                         NA,
                                         body.variance.full$mean.female)

body.variance.full$n.female <- ifelse(body.variance.full$n.female==1,
                                      NA,
                                      body.variance.full$n.female)

body.variance.full$SD.female <- ifelse(body.variance.full$n.female==1,
                                       NA,
                                       body.variance.full$SD.female)


# assigning NA to mean and n to any estimate with a missing SD
body.variance.full$mean.male <- ifelse(is.na(body.variance.full$SD.male),
                                       NA,
                                       body.variance.full$mean.male)

body.variance.full$n.male <- ifelse(is.na(body.variance.full$SD.male),
                                    NA,
                                    body.variance.full$n.male)

body.variance.full$mean.female <- ifelse(is.na(body.variance.full$SD.female),
                                         NA,
                                         body.variance.full$mean.female)

body.variance.full$n.female <- ifelse(is.na(body.variance.full$SD.female),
                                      NA,
                                      body.variance.full$n.female)


# assigning NA to mean, SD and n to any estimate with SD == 0
body.variance.full$mean.male <- ifelse(body.variance.full$SD.male==0,
                                       NA,
                                       body.variance.full$mean.male)

body.variance.full$n.male <- ifelse(body.variance.full$SD.male==0,
                                    NA,
                                    body.variance.full$n.male)

body.variance.full$SD.male <- ifelse(body.variance.full$SD.male==0,
                                     NA,
                                     body.variance.full$SD.male)

body.variance.full$mean.female <- ifelse(body.variance.full$SD.female==0,
                                         NA,
                                         body.variance.full$mean.female)

body.variance.full$n.female <- ifelse(body.variance.full$SD.female==0,
                                      NA,
                                      body.variance.full$n.female)

body.variance.full$SD.female <- ifelse(body.variance.full$SD.female==0,
                                       NA,
                                       body.variance.full$SD.female)


################################################################################
# plotting mean wing length female ~ male with cleaned dataset
#par(mfrow = c(1, 2))

# plot(body.variance.full$mean.male,body.variance.full$mean.female,
#      main = "Wing length Female~Male Overall correlation",
#      xlab = "Mean Wing length (Males)",
#      ylab = "Mean Wing length (Females)",
#      xlim = c(0,605),
#      ylim = c(0,605))
# abline(1,1)

# and estimating the correlation
cor.test(body.variance.full$mean.male,body.variance.full$mean.female)


################################################################################
# Figures for: S1. Male-Female wing length correlations
################################################################################

# https://github.com/IndrajeetPatil/ggstatsplot/blob/main/README.md
p1 <- ggscatterstats(
  data = body.variance.full, ## data frame from which variables are taken
  x = mean.male, ## predictor/independent variable
  y = mean.female, ## dependent variable
  type = "parametric",
  digits = 2L,
  bf.message = F, # remove Bayesian test
  point.args = list(alpha = 0.25,size = 3,stroke = 1,color="blue4"),
  point.width.jitter = 0,
  point.height.jitter = 0,
  smooth.line.args = list(linewidth = 1.25, 
                          color = "black", 
                          method = "lm", 
                          formula = y ~ x),
  xsidehistogram.args = list(fill = "#009E73", color = "black", na.rm = TRUE,
                             binwidth = 10),
  ysidehistogram.args = list(fill = "#D55E00", color = "black", na.rm = TRUE,
                             binwidth = 10),  xlab = "Mean Wing length [Males]", ## label for the x-axis
  ylab = "Mean Wing length [Females]", ## label for the y-axis
  label.var = species, ## variable to use for labeling data points
  label.expression = mean.male > 550 & mean.female < 450, ## expression for deciding which points to label
  point.label.args = list(alpha = 0.7, size = 4, color = "grey25"),
  # xfill = "#CC79A7", ## fill for marginals on the x-axis
  # yfill = "#009E73", ## fill for marginals on the y-axis
  #title = "Wing length Female ~ Male Overall correlation",
  #caption = ""
) + ## modifying the plot further
  ggplot2::scale_y_continuous(
    breaks = seq(from = 0, to = 600, by = 50)
  )

p1

# then doing the same by using log-transformed values
# plot(log(body.variance.full$mean.male),log(body.variance.full$mean.female),
#      main = "Wing length Female~Male Overall correlation",
#      xlab = "ln(Mean Wing length (Males))",
#      ylab = "ln(Mean Wing length (Females))",
#      xlim = c(3.5,6.5),
#      ylim = c(3.5,6.5))
# abline(0,1)

# and estimating the correlation
cor.test(log(body.variance.full$mean.male),log(body.variance.full$mean.female))

body.variance.full$mean.male.log <- log(body.variance.full$mean.male)
body.variance.full$mean.female.log <- log(body.variance.full$mean.female)

# https://github.com/IndrajeetPatil/ggstatsplot/blob/main/README.md
p2 <- ggscatterstats(
  data = body.variance.full, ## data frame from which variables are taken
  x = mean.male.log, ## predictor/independent variable
  y = mean.female.log, ## dependent variable
  type = "parametric",
  digits = 2L,
  bf.message = F, # remove Bayesian test
  point.args = list(alpha = 0.25,size = 3,stroke = 1,color="blue4"),
  point.width.jitter = 0,
  point.height.jitter = 0,
  smooth.line.args = list(linewidth = 1.25, 
                          color = "black", 
                          method = "lm", 
                          formula = y ~ x),
  xsidehistogram.args = list(fill = "#009E73", color = "black", na.rm = TRUE,
                             binwidth = 0.05),
  ysidehistogram.args = list(fill = "#D55E00", color = "black", na.rm = TRUE,
                             binwidth = 0.05),
  xlab = "ln (Mean Wing length [Males])", ## label for the x-axis
  ylab = "ln (Mean Wing length [Females])", ## label for the y-axis
  label.var = species, ## variable to use for labeling data points
  label.expression = mean.male > 550 & mean.female < 450, ## expression for deciding which points to label
  point.label.args = list(alpha =0.5, size = 4, color = "grey25"),
  # xfill = "#CC79A7", ## fill for marginals on the x-axis
  # yfill = "#009E73", ## fill for marginals on the y-axis
  #title = "Wing length Female ~ Male Overall correlation",
  #caption = ""
) + ## modifying the plot further
  ggplot2::scale_y_continuous(
    breaks = seq(from = 3.5, to = 6.5, by = 0.2)
  )

p2

## combining the individual plots into a single plot: Figure S1
combine_plots(
  list(p1, p2),
  plotgrid.args = list(nrow = 1),
  annotation.args = list(
    title = "Wing length Female ~ Male Overall correlation",
    caption = ""
  )
)


################################################################################
# Do we have evidence for Taylor's law?
# S2. Mean-Variance relationship and CV summary statistics
################################################################################

#par(mfrow = c(1, 2))

# males
# plot(log(body.variance.full$mean.male),log(body.variance.full$SD.male),
#      main = "Males: Wing length mean~SD Overall correlation",
#      xlab = "ln(mean)",
#      ylab = "ln(SD)",
#      xlim = c(3.5,6.5),
#      ylim = c(-1,4))
# 
# abline(-3.75,1)
# 
# regression to add line to the plot
taylors.male <- lm(log(body.variance.full$SD.male)~log(body.variance.full$mean.male))
summary(taylors.male)
confint(taylors.male)
# 
# abline(taylors.male$coefficients[1],taylors.male$coefficients[2],col = "red")
# 
# # and estimating the correlation
# cor.test(log(body.variance.full$mean.male),log(body.variance.full$SD.male))

body.variance.full$SD.male.log <- log(body.variance.full$SD.male)

# https://github.com/IndrajeetPatil/ggstatsplot/blob/main/README.md
# https://indrajeetpatil.github.io/ggstatsplot/reference/ggscatterstats.html
p3 <- ggscatterstats(
  data = body.variance.full, ## data frame from which variables are taken
  x = mean.male.log, ## predictor/independent variable
  y = SD.male.log, ## dependent variable
  type = "parametric",
  digits = 2L,
  bf.message = F, # remove Bayesian test
  point.args = list(alpha = 0.25,
                    #size = 3,
                    #size = body.variance.full$n.male,
                    stroke = 1, 
                    color="blue4",
                    aes(size=n.male)),
  point.width.jitter = 0,
  point.height.jitter = 0,
  smooth.line.args = list(linewidth = 1.25, 
                          color = "black", 
                          method = "lm", 
                          formula = y ~ x),
  xsidehistogram.args = list(fill = "#009E73", color = "black", na.rm = TRUE,
                             binwidth = 0.05),
  ysidehistogram.args = list(fill = "#D55E00", color = "black", na.rm = TRUE,
                             binwidth = 0.1),  
  xlab = "ln (mean) [Male]", ## label for the x-axis
  ylab = "ln (SD) [Male]", ## label for the y-axis
  label.var = species, ## variable to use for labeling data points
  #label.expression = mean.male.log < 4.1 & SD.male.log < -0.15, ## expression for deciding which points to label
  label.expression = n.male < 4, ## expression for deciding which points to label
  point.label.args = list(alpha = 0.7, size = 4, color = "grey25"),
) + ## modifying the plot further
  ggplot2::scale_y_continuous(
    breaks = seq(from = -1.1, to = 3.7, by = 0.3)
  ) + ## modifying the plot further
  ggplot2::scale_x_continuous(
    breaks = seq(from = 3.5, to = 6.5, by = 0.3)
  ) + ggplot2::theme(legend.position="none")

p3


# females
# plot(log(body.variance.full$mean.female),log(body.variance.full$SD.female),
#      main = "Females: Wing length mean~SD Overall correlation",
#      xlab = "ln(mean)",
#      ylab = "ln(SD)",
#      xlim = c(3.5,6.5),
#      ylim = c(-1,4))
# 
# abline(-3.75,1)

# regression to add line to the plot
taylors.female <- lm(log(body.variance.full$SD.female)~log(body.variance.full$mean.female))
summary(taylors.female)
confint(taylors.female)

#abline(taylors.female$coefficients[1],taylors.female$coefficients[2],col = "red")

# and estimating the correlation
cor.test(log(body.variance.full$mean.female),log(body.variance.full$SD.female))


body.variance.full$SD.female.log <- log(body.variance.full$SD.female)

# https://github.com/IndrajeetPatil/ggstatsplot/blob/main/README.md
# https://indrajeetpatil.github.io/ggstatsplot/reference/ggscatterstats.html
p4 <- ggscatterstats(
  data = body.variance.full, ## data frame from which variables are taken
  x = mean.female.log, ## predictor/independent variable
  y = SD.female.log, ## dependent variable
  type = "parametric",
  digits = 2L,
  bf.message = F, # remove Bayesian test
  point.args = list(alpha = 0.25,
                    #size = 3,
                    #size = body.variance.full$n.male,
                    stroke = 1, 
                    color="blue4",
                    aes(size=n.female)),
  point.width.jitter = 0,
  point.height.jitter = 0,
  smooth.line.args = list(linewidth = 1.25, 
                          color = "black", 
                          method = "lm", 
                          formula = y ~ x),
  xsidehistogram.args = list(fill = "#009E73", color = "black", na.rm = TRUE,
                             binwidth = 0.05),
  ysidehistogram.args = list(fill = "#D55E00", color = "black", na.rm = TRUE,
                             binwidth = 0.1),  
  xlab = "ln (mean) [Female]", ## label for the x-axis
  ylab = "ln (SD) [Female]", ## label for the y-axis
  label.var = species, ## variable to use for labeling data points
  #label.expression = mean.male.log < 4.1 & SD.male.log < -0.15, ## expression for deciding which points to label
  label.expression = n.female < 4, ## expression for deciding which points to label
  point.label.args = list(alpha = 0.7, size = 4, color = "grey25"),
) + ## modifying the plot further
  ggplot2::scale_y_continuous(
    breaks = seq(from = -1.1, to = 3.7, by = 0.3)
  ) + ## modifying the plot further
  ggplot2::scale_x_continuous(
    breaks = seq(from = 3.5, to = 6.5, by = 0.3)
  ) + ggplot2::theme(legend.position="none")

p4


## combining the individual plots into a single plot: Figure S2
combine_plots(
  list(p3, p4),
  plotgrid.args = list(nrow = 1),
  annotation.args = list(
    title = "Taylor's Law",
    caption = ""
  )
)


################################################################################
# reference
table(body.variance.full$reference)


################################################################################
# re-organizing dataset for the analyses: males vs female
################################################################################
# sorting database alphabetically
body.variance.full <- body.variance.full[order(body.variance.full$species),]


# creating an identifier for those male and female estimates that were at least 
# reported together, since we are going to use female and male data as 
# independent data points from now on
body.variance.full$pair_ID <- 1:nrow(body.variance.full)


# splitting the database so that males and females can be analyzed together
body.variance.full.final.males <- as.data.frame(body.variance.full 
                                                # removing female estimates
                                                %>% select(-mean.female, 
                                                           -n.female, 
                                                           -SD.female) 
                                                # excluding NA values (including 
                                                # those assigned above)
                                                %>% filter(!is.na(mean.male))
                                                # renaming variables
                                                %>% rename(mean = mean.male,
                                                           n = n.male, 
                                                           SD = SD.male)
                                                # categorizing estimates by sex
                                                %>% mutate(sex = "male"))

nrow(body.variance.full.final.males)


body.variance.full.final.females <- as.data.frame(body.variance.full 
                                                  # removing male estimates
                                                  %>% select(-mean.male, 
                                                             -n.male, 
                                                             -SD.male)
                                                  # excluding NA values (including 
                                                  # those assigned above)
                                                  %>% filter(!is.na(mean.female))
                                                  # renaming variables
                                                  %>% rename(mean = mean.female,
                                                             n = n.female,
                                                             SD = SD.female)
                                                  # categorizing estimates by sex
                                                  %>% mutate(sex = "female"))


nrow(body.variance.full.final.females)


# building the combined male-female database for the analyses
body.variance.full.sex.sep <- rbind(body.variance.full.final.males,
                                    body.variance.full.final.females)

# making sex a factor
body.variance.full.sex.sep$sex <- as.factor(body.variance.full.sex.sep$sex)

# exploring the new combined dataset
summary(body.variance.full.sex.sep)

# How many studies in total?
length(unique(body.variance.full.sex.sep$source))
length(unique(body.variance.full.sex.sep$reference_link))

# How many populations in total?
length(unique(body.variance.full.sex.sep$population_ID))


################################################################################
# Phylogeny reconstruction and species list fixing
################################################################################

# obtaining dataframe listing the Open Tree identifiers potentially matching
# our list of species (be aware that this will take a few minutes. You can avoid
# waiting by loading the data below)
study.species.original <- as.character(unique(body.variance.full.sex.sep$species))

# How many (sub)species?
length(unique(study.species.original))

# # extracting taxonomic information
# taxa <- tnrs_match_names(names = study.species.original)
# 
# # check approximate matches
# taxa[taxa$approximate_match==TRUE & !(is.na(taxa$approximate_match)),]
# 
# # check approximate matches
# taxa[taxa$is_synonym==TRUE & !(is.na(taxa$is_synonym)),]
# 
# # check number of matches
# taxa[taxa$number_matches>1 & !(is.na(taxa$number_matches)),]

# Based on the checks above, we had to implement some changes to the species names:
# The remaining changes are typos that needed correction, update in species names or bringing subspecies up to species level so that everything is comparable
body.variance.full.sex.sep$species.updated <-
  recode(body.variance.full.sex.sep$species,
         "Ammodramus savannarum ammolegus" = "Ammodramus savannarum", #not matched at the tree of life so brought down to species level
         "Ammodramus savannarum perpallidus" = "Ammodramus savannarum", #not matched at the tree of life so brought down to species level
         "Ammodramus savannarum pratensis" = "Ammodramus savannarum", #not matched at the tree of life so brought down to species level
         "Calidris alpina schinzii" = "Calidris alpina", #not matched at the tree of life so brought down to species level
         "Carpodacus mexicanus frontalis" = "Haemorhous mexicanus", #not matched at the tree of life so brought down to species level
         "Malurus alboscapulatus alboscapulatus" = "Malurus alboscapulatus", #not matched at the tree of life so brought down to species level
         "Malurus alboscapulatus kutubu" = "Malurus alboscapulatus", #not matched at the tree of life so brought down to species level
         "Malurus alboscapulatus naimii (highland form)" = "Malurus alboscapulatus", #not matched at the tree of life so brought down to species level
         "Malurus alboscapulatus naimii (lowland form)" = "Malurus alboscapulatus", #not matched at the tree of life so brought down to species level
         "Malurus melanocephalus cruentatus" = "Malurus melanocephalus", #not matched at the tree of life so brought down to species level
         "Miliaria calandra" = "Emberiza calandra", #not matched at the tree of life so name updated
         
         "Accipter cooperii" = "Accipiter cooperii", #typo fixed
         "Accipter cooperi" = "Accipiter cooperii", #typo fixed
         "Larus crassirostis" = "Larus crassirostris", #typo fixed
         "Marmoronetta angustirostris" = "Marmaronetta angustirostris", #typo fixed
         "Onychoprion fuscata" = "Onychoprion fuscatus", #typo fixed
         "Phaethon rubricaudus" = "Phaethon rubricauda", #typo fixed
         #"Phaeton aethereus" = "Phaethon aethereus", #typo fixed
         #"Phaeton lepturus" = "Phaethon lepturus", #typo fixed
         
         # subspecies to species
         "Porphyrio porphyrio" = "Porphyrio porphyrio porphyrio",
         
         "Andropadus latirostris" = "Eurillas latirostris", #updating name
         "Andropadus tephrolaemus" = "Arizelocichla tephrolaema", #updating name
         "Bonasa sewerzowi" = "Tetrastes sewerzowi", #updating name
         "Carduelis cannabina" = "Linaria cannabina", #updating name
         "Carduelis chloris" = "Chloris chloris", #updating name
         "Carpodacus mexicanus" = "Haemorhous mexicanus", #updating name
         "Catharacta skua" = "Stercorarius skua", #updating name
         "Chlidonias hybrida" = "Chlidonias hybridus", #updating name
         "Dendrocopos pubescens" = "Picoides pubescens", #updating name
         "Diglossa cyanea" = "Diglossopis cyanea", #updating name
         "Larus ridibundus" = "Chroicocephalus ridibundus", #updating name
         "Myiothlypis luteoviridis" = "Basileuterus luteoviridis", #updating name
         "Myiothlypis signata" = "Basileuterus signatus", #updating name
         "Parus caeruleus" = "Cyanistes caeruleus", #updating name
         "Parus caeruleus ogliastrae" = "Cyanistes caeruleus", #updating name and bringing up to species to make everything comparable
         "Philomachus pugnax" = "Calidris pugnax", #updating name
         "Phylloscopus inornatus" = "Abrornis inornata", #updating name
         "Spilopelia senegalensis" = "Streptopelia senegalensis", #updating name
         "Cyanopica cyanus" = "Cyanopica cooki", #updating name based on population of origin
         
         "Serinus citrinella" = "Carduelis citrinella", #changed the name because otherwise Serinus citrinella seemingly assigns it to what is known now as a different species living in Corsica
         
         # more subsp to sp
         "Ammodramus savannarum floridanus" = "Ammodramus savannarum", #brought down to species level because the other three subsp were not matched in the open tree of life. This way, we make the floridanus subsp comparable to the other three subsp
         "Malurus alboscapulatus aida" = "Malurus alboscapulatus", #brought down to species level because the other subsp were not matched in the open tree of life. This way, we make this subsp comparable to the other subsp in the genus
         "Malurus alboscapulatus balim" = "Malurus alboscapulatus", #brought down to species level because the other subsp were not matched in the open tree of life. This way, we make this subsp comparable to the other subsp in the genus
         "Malurus alboscapulatus lorentzi" = "Malurus alboscapulatus", #brought down to species level because the other subsp were not matched in the open tree of life. This way, we make this subsp comparable to the other subsp in the genus
         "Malurus alboscapulatus randi" = "Malurus alboscapulatus", #brought down to species level because the other subsp were not matched in the open tree of life. This way, we make this subsp comparable to the other subsp in the genus
         "Malurus alboscapulatus moretoni" = "Malurus alboscapulatus", #brought down to species level because the other subsp were not matched in the open tree of life. This way, we make this subsp comparable to the other subsp in the genus
         "Malurus melanocephalus melanocephalus" = "Malurus melanocephalus", #brought down to species level because the other subsp were not matched in the open tree of life. This way, we make this subsp comparable to the other subsp in the genus
         "Acanthis flammea islandica" = "Acanthis flammea", # although it could probably stay as subspecies and nothing would change
         #"Anas cyanoptera septentrionalium" = "Anas cyanoptera", # not brought up to species because this species is no longer monophyletic, which means that only by keeping it as subspecies we can get it's phylogenetic position
         "Calonectris diomedea diomedea" = "Calonectris diomedea", # although it could probably stay as subspecies and nothing would change
         "Emberiza schoeniclus lusitanica" = "Emberiza schoeniclus", # although it could probably stay as subspecies and nothing would change
         "Emberiza schoeniclus schoeniclus" = "Emberiza schoeniclus", # although it could probably stay as subspecies and nothing would change
         "Emberiza schoeniclus witherbyi" = "Emberiza schoeniclus", # although it could probably stay as subspecies and nothing would change
         "Haematopus ostralegus longipes" = "Haematopus ostralegus", # although it could probably stay as subspecies and nothing would change
         "Haematopus ostralegus ostralegus" = "Haematopus ostralegus", # although it could probably stay as subspecies and nothing would change
         "Hydrobates pelagicus pelagicus" = "Hydrobates pelagicus", #brought down to species level because there are other instances with Hydrobates pelagicus only, just to make it comparable
         "Malurus leucopterus edouardi" = "Malurus leucopterus", # brought down to species because subsp ranges largely overlap and also because this was the only species with multiple subspecies left, meaning that bringing it down to the species level makes it more comparable to the rest of the dataset
         "Malurus leucopterus leuconotus" = "Malurus leucopterus", # brought down to species because subsp ranges largely overlap and also because this was the only species with multiple subspecies left, meaning that bringing it down to the species level makes it more comparable to the rest of the dataset
         "Malurus leucopterus leucopterus" = "Malurus leucopterus", # brought down to species because subsp ranges largely overlap and also because this was the only species with multiple subspecies left, meaning that bringing it down to the species level makes it more comparable to the rest of the dataset
         "Zonotrichia leucophrys oriantha" = "Zonotrichia leucophrys", # although it could probably stay as subspecies and nothing would change
         "Gallinula galeata sandvicensis" = "Gallinula galeata",
         .default = levels(body.variance.full.sex.sep$species))

# making sepcies.updated a factor
body.variance.full.sex.sep$species.updated <- factor(body.variance.full.sex.sep$species.updated)

# How many (sub)species?
length(unique(body.variance.full.sex.sep$species.updated))


# summary statistics:
# males
sort(table(body.variance.full.sex.sep[body.variance.full.sex.sep$sex=="male",
                                      "species.updated"]))
length(unique(body.variance.full.sex.sep[body.variance.full.sex.sep$sex=="male",
                                         "species.updated"]))
summary(as.numeric(table(body.variance.full.sex.sep[body.variance.full.sex.sep$sex=="male",
                                                    "species.updated"])))
sd(as.numeric(table(body.variance.full.sex.sep[body.variance.full.sex.sep$sex=="male",
                                               "species.updated"])))

# sample sizes
summary(body.variance.full.sex.sep[body.variance.full.sex.sep$sex=="male",
                                   c("species.updated","n")])


# summary statistics:
# females
sort(table(body.variance.full.sex.sep[body.variance.full.sex.sep$sex=="female",
                                      "species.updated"]))
length(unique(body.variance.full.sex.sep[body.variance.full.sex.sep$sex=="female",
                                         "species.updated"]))
# mean
summary(
  as.numeric(
    table(
      body.variance.full.sex.sep[body.variance.full.sex.sep$sex=="female","species.updated"]
    )
  )[as.numeric(table(body.variance.full.sex.sep[body.variance.full.sex.sep$sex=="female",
                                                "species.updated"]))!=0]) # excluding 0 assigned to species for which there is no female data

sd(
  as.numeric(
    table(
      body.variance.full.sex.sep[body.variance.full.sex.sep$sex=="female","species.updated"]
    )
  )[as.numeric(table(body.variance.full.sex.sep[body.variance.full.sex.sep$sex=="female",
                                                "species.updated"]))!=0]) # excluding 0 assigned to species for which there is no female data

# sample sizes
summary(body.variance.full.sex.sep[body.variance.full.sex.sep$sex=="female",
                                   c("species.updated","n")])


# repeating the taxonomic extraction process after fixing the species names
study.species <- as.character(unique(body.variance.full.sex.sep$species.updated))
# taxa <- tnrs_match_names(names = study.species) # no warnings this time, good sign
# taxa[taxa$approximate_match==TRUE & !(is.na(taxa$approximate_match)),] # none
# taxa[taxa$is_synonym==TRUE & !(is.na(taxa$is_synonym)),] # none
# taxa[taxa$number_matches>1 & !(is.na(taxa$number_matches)),] # none
# # all fixed now
# taxa
# 
# # some further checks: no longer needed
# ott_id_tocheck <- taxa[taxa$number_matches != 1,"ott_id"]
# ott_id_tocheck
# 
# # no longer needed because no species with number_matches > 1. Nice.
# for(i in 1:length(ott_id_tocheck)){
#   print(inspect(taxa, ott_id = ott_id_tocheck[i]))
# }
# 
# # saving the taxonomic data last created on ... to speed the process in the
# # future and allow full computational reproducibility
# save(taxa,file = "data/phylogeny/taxa_Open_Tree_of_Life_20240909.RData")

# loading the most updated taxonomic data created
load("data/phylogeny/taxa_Open_Tree_of_Life_20240909.RData") #taxa
# everything seems in order now

# # retrieving phylogenetic relationships among taxa in the form of a trimmed
# # sub-tree
# tree <- tol_induced_subtree(ott_ids = taxa[["ott_id"]], label_format = "name")
# 
# # we need to check for the existence of polytomies
# is.binary(tree)
# # there used to be a polytomy for Anous stolidus and Chroicocephalus ridibundus
# # which no longer exist
# # (see https://doi.org/10.1111/j.1469-7998.2009.00622.x)
# # but there is one polytomy that exist for Nectarina famosa and Nectarinia johnstoni
# 
# # to take care of this polytomy, we are going to use a randomization approach
# set.seed(23)
# tree_random <- multi2di(tree,random=TRUE)
# is.binary(tree_random)
# 
# # to confirm that our tree covers all the species we wanted to include, and
# # make sure that the species names in our database match those in the tree, we
# # use the following code
# 
# # here are the species included in the tree
# sort(tree_random$tip.label)
# 
# # removing the underscore "_" from the tree tip.label
# tree_random$tip.label <- gsub("_"," ", tree_random$tip.label)
# sort(tree_random$tip.label)
# 
# # here are the species included in our dataset
# sort(unique(as.character(body.variance.full.sex.sep$species.updated)))
# 
# # comparing the tree to our dataset
# intersect(as.character(tree_random$tip.label),
#           as.character(body.variance.full.sex.sep$species.updated))
# setdiff(as.character(body.variance.full.sex.sep$species.updated),
#         as.character(tree_random$tip.label)) # listed in our database but not in the tree
# setdiff(as.character(tree_random$tip.label),
#         as.character(body.variance.full.sex.sep$species.updated)) # listed in the tree but not in our database
# 
# # we can now save the tree
# save(tree_random, file = "data/phylogeny/tree_random_20240909.Rdata")
# 
# 
# # compute branch lengths of tree
# phylo_branch <- compute.brlen(tree_random, method = "Grafen", power = 1)
# 
# # check tree is ultrametric
# is.ultrametric(phylo_branch) # TRUE
# 
# # matrix to be included in the models
# phylo_cor <- vcv(phylo_branch, cor = T)
# 
# # finally, save matrix for future analyses
# save(phylo_cor, file = "data/phylogeny/phylo_cor_20240909.Rdata")


# loading tree and matrix
load("data/phylogeny/tree_random_20240909.Rdata")
load("data/phylogeny/phylo_cor_20240909.Rdata")


# # we can then plot the tree
# plot(tree_random, type = "fan", cex=0.8, label.offset =.1, no.margin = TRUE) #check: https://www.rdocumentation.org/packages/ape/versions/5.3/topics/plot.phylo
# 
# plot(tree_random, cex=0.25, label.offset =.1, no.margin = TRUE) #check: https://www.rdocumentation.org/packages/ape/versions/5.3/topics/plot.phylo


################################################################################
# Preparing template to obtain migratory and feeding type information
################################################################################

# first creating a database to be exported and filled in with data
species.database <- unique(body.variance.full.sex.sep[,c("species.updated",
                                                         "location")])

# adding variables needed to be filled up
species.database$migratory <- ""
species.database$migratory.comments <- ""
species.database$feeding.type <- ""
species.database$feeding.type.comments <- ""
species.database$broader.migration <- ""
species.database$broader.feeding.type <- ""
species.database$broader.comments <- ""

species.database <- unique(species.database[,c("species.updated",
                                               "migratory",
                                               "migratory.comments",
                                               "feeding.type",
                                               "feeding.type.comments",
                                               "broader.migration",
                                               "broader.feeding.type",
                                               "broader.comments",
                                               "location")])

species.database <- species.database[order(species.database$species.updated),]

# # saving database so that AST and mostly KR can fill this information up for the 
# # analyses
# write.csv(unique(species.database[,c("species.updated",
#                                      "migratory",
#                                      "migratory.comments",
#                                      "feeding.type",
#                                      "feeding.type.comments",
#                                      "broader.migration",
#                                      "broader.feeding.type",
#                                      "broader.comments")]),
#           "data/final/99_final_bird_species_list_extraction_template.csv",
#           row.names=FALSE)


################################################################################
# Adding the species-specific and population-specific information about 
# migration and feeding type
################################################################################

# importing data
species.information <- read.csv("data/final/02_final_bird_species_list_migration_and_feeding_vagrant.csv",
                                header=T, sep=",")

# but before, we need to explore the dataset to make sure we do not miss any 
# information and/or that the information we have is consistent
species.information[order(species.information$species.updated),
                    c("species.updated", "migratory", "feeding.type",
                      "broader.migration","broader.feeding.type")]

# reformatting some variables types to factors
cols.factors.2 <- c("species.updated", "migratory", "feeding.type",
                    "broader.migration","broader.feeding.type")
species.information[cols.factors.2] <- lapply(species.information[cols.factors.2], 
                                              factor)


################################################################################
# checking for potential missing species

# species names in species.information but not in body.variance.full.sex.sep
setdiff(species.information$species.updated,
        unique(body.variance.full.sex.sep$species.updated))

# species names in body.variance.full.sex.sep but not in species.information
setdiff(unique(body.variance.full.sex.sep$species.updated),
        species.information$species.updated)

# exploring species for which we missed information
missing.species <- body.variance.full.sex.sep[body.variance.full.sex.sep$species.updated %in%
                                                setdiff(unique(body.variance.full.sex.sep$species.updated),
                                                        species.information$species.updated),
                                              c("species","species.updated",
                                                "source","reference_link",
                                                "location","sex")]
missing.species <- missing.species[order(missing.species$species.updated),]
missing.species

# exporting the list to search for the missing information
# write.csv(missing.species[missing.species$sex=="male",],
#           "data/missing_species_information.csv",
#           row.names=FALSE)


################################################################################
# searching for potential disagreements between feeding.type and
# broader.feeding.type
species.information[species.information$feeding.type=="no" & 
                      species.information$broader.feeding.type=="yes",]

species.information[species.information$feeding.type=="yes" & 
                      species.information$broader.feeding.type=="no",]

# searching for potential disagreements between migratory and broader.migration
species.information[species.information$migratory=="no" & 
                      species.information$broader.migration=="yes",]

species.information[species.information$migratory=="no" & 
                      species.information$broader.migration=="partially",]

species.information[species.information$migratory=="yes" & 
                      species.information$broader.migration=="no",]


################################################################################
# Finally, merging databases to create the final database for analyses
body.variance.full.sex.sep.species.info <- merge(body.variance.full.sex.sep,
                                                 species.information,
                                                 all.x = T, 
                                                 by = "species.updated")


# # last fix to the dataset, removing an estimate based on three individuals 
# # which has an SD==0
# body.variance.full.sex.sep.species.info <- 
#   as.data.frame(body.variance.full.sex.sep.species.info %>% filter(SD > 0))


################################################################################
# adding population-specific migratory status for those species showing both
# migration yes and no

# importing population-specific cases
population.specific.cases <- read.csv("data/final/03_population-specific_migratory_status.csv",
                                      header=T, sep=",")
population.specific.cases


# how does entries look like
body.variance.full.sex.sep.species.info[body.variance.full.sex.sep.species.info$species
                                        %in% (population.specific.cases$species.updated),
                                        c("species.updated","source","population_ID",
                                          "location")]

# modifying the current "partially" value in "broader.migration" using the 
# information contained in the population-specific cases database

# first adding the value to the current dataset
body.variance.full.sex.sep.species.info <- merge(body.variance.full.sex.sep.species.info,
                                                 unique(population.specific.cases[,c("population_ID","broader.migration.2")]),
                                                 all.x = T)

# then substituting "partially" values
table(body.variance.full.sex.sep.species.info$broader.migration)
table(is.na(body.variance.full.sex.sep.species.info$broader.migration))
body.variance.full.sex.sep.species.info$broader.migration <- as.character(body.variance.full.sex.sep.species.info$broader.migration)
body.variance.full.sex.sep.species.info$broader.migration.2 <- as.character(body.variance.full.sex.sep.species.info$broader.migration.2)

body.variance.full.sex.sep.species.info$broader.migration <- ifelse(body.variance.full.sex.sep.species.info$broader.migration=="partially",
                                                                    body.variance.full.sex.sep.species.info$broader.migration.2,
                                                                    body.variance.full.sex.sep.species.info$broader.migration)

table(body.variance.full.sex.sep.species.info$broader.migration)
table(is.na(body.variance.full.sex.sep.species.info$broader.migration))

# removing population-specific variable
body.variance.full.sex.sep.species.info$broader.migration.2 <- NULL
names(body.variance.full.sex.sep.species.info)
body.variance.full.sex.sep.species.info$broader.migration <- as.factor(body.variance.full.sex.sep.species.info$broader.migration)
summary(body.variance.full.sex.sep.species.info)


################################################################################
# Dealing with those classified as vagrants

# we are categorising them as "non-migratory" for the main analysis, but 
# we will consider them as "migrants" for a sensitivity analysis

# creating two new variables where vagrants are considered migrants
body.variance.full.sex.sep.species.info$broader.migration.inc.vagrants <- ifelse(as.character(body.variance.full.sex.sep.species.info$broader.migration)=="vagrant",
                                                                                 "yes",
                                                                                 as.character(body.variance.full.sex.sep.species.info$broader.migration))

body.variance.full.sex.sep.species.info$migratory.inc.vagrants <- ifelse(as.character(body.variance.full.sex.sep.species.info$broader.migration)=="vagrant",
                                                                         "yes",
                                                                         as.character(body.variance.full.sex.sep.species.info$migratory))

# then substituting "vagrants" by "no" in broader.migration, which is the
# variable that will be used for the main analyses
body.variance.full.sex.sep.species.info$broader.migration <- ifelse(as.character(body.variance.full.sex.sep.species.info$broader.migration)=="vagrant",
                                                                    "no",
                                                                    as.character(body.variance.full.sex.sep.species.info$broader.migration))


################################################################################
# Exporting final clean dataset for analyses
################################################################################

# renaming population level columns before saving
body.variance.full.sex.sep.species.info <- rename(body.variance.full.sex.sep.species.info, 
                                                  migratory.pop.level = broader.migration)

body.variance.full.sex.sep.species.info <- rename(body.variance.full.sex.sep.species.info, 
                                                  migratory.inc.vagrants.pop.level = broader.migration.inc.vagrants)

body.variance.full.sex.sep.species.info <- rename(body.variance.full.sex.sep.species.info,
                                                  feeding.type.pop.level = broader.feeding.type)

body.variance.full.sex.sep.species.info <- rename(body.variance.full.sex.sep.species.info,
                                                  study_ID = source)

body.variance.full.sex.sep.species.info <- rename(body.variance.full.sex.sep.species.info,
                                                  subset = reference)

body.variance.full.sex.sep.species.info <- rename(body.variance.full.sex.sep.species.info,
                                                  measurement_notes = Measurement)


# removing unnecessary variables
body.variance.full.sex.sep.species.info$exclude <- NULL

# data based used in further scripts
write.csv(body.variance.full.sex.sep.species.info,
          "data/final/04_final_full_and_clean_bird_wing_size_dataset.csv",
          row.names=FALSE)
