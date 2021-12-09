# Author's information: 
# Alfredo Sanchez-Tojar (@ASanchez_Tojar)
# Profile: https://goo.gl/PmpPEB
# Department of Evolutionary Biology, Bielefeld University (GER) 
# Email: alfredo.tojar@gmail.com
# Script first created on the 25th of August 2021

pacman::p_load(ape,fulltext,metafor,rotl, #openxlsx,readxl,
               treebase,diagram,dplyr,stringr,ggplot2,
               tidyverse,ggpubr,orchaRd) #pander,brms,

# cleaning
rm(list=ls())

#function
### Functions for I2
I2 <- function(model, method = c("Wolfgang", "Shinichi")){
  warning("Make sure you have the observation (effec size) level random effect\n")
  ## evaluate choices
  method <- match.arg(method)
  
  # Wolfgang's method
  if(method == "Wolfgang"){
    W <- solve(model$V) 
    X <- model.matrix(model)
    P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    I2_total <- sum(model$sigma2) / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P)))
    I2_each  <- model$sigma2 / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P)))
    names(I2_each) = paste0("I2_", model$s.names)
    
    # putting all together
    I2s <- c(I2_total = I2_total, I2_each)
    
    # or my way
  } else {
    # sigma2_v = typical sampling error variance
    sigma2_v <- sum(1/model$vi) * (model$k-1) / (sum(1/model$vi)^2 - sum((1/model$vi)^2)) 
    I2_total <- sum(model$sigma2) / (sum(model$sigma2) + sigma2_v) #s^2_t = total variance
    I2_each  <- model$sigma2 / (sum(model$sigma2) + sigma2_v)
    names(I2_each) = paste0("I2_", model$s.names)
    
    # putting all together
    I2s <- c(I2_total = I2_total, I2_each)
  }
  return(I2s)
}

### Functions for R2

R2 <- function(model){
  warning("Make sure you have the observation (effec size) level random effect as the last in the formula\n")
  
  # fixed effect variance
  fix <- var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))
  
  # marginal
  R2m <- fix / (fix + sum(model$sigma2))
  R2
  #Rm <- round(100*R2m, 3)
  
  # conditional
  R2c <- (fix + sum(model$sigma2) - model$sigma2[length(model$sigma2)]) / 
    (fix + sum(model$sigma2))
  
  R2s <- c(R2_marginal = R2m, R2_coditional = R2c)
  
  return(R2s)
}


# extracting mean and CI from each metafor model
estimates.CI <- function(model){
  db.mf <- data.frame(model$b,row.names = 1:nrow(model$b))
  db.mf <- cbind(db.mf,model$ci.lb,model$ci.ub,row.names(model$b))
  names(db.mf) <- c("mean","lower","upper","estimate")
  return(db.mf[,c("estimate","mean","lower","upper")])
}


# importing data
body.variance.full <- read.csv("data/final/final_bird_wing_size_dataset.csv",header=T,sep=",")

# reformatting some variables types to factors
cols.factors <- c("species","source","reference_link","population_ID","reference","trait","unit","exclude","location","Measurement")
body.variance.full[cols.factors] <- lapply(body.variance.full[cols.factors], factor)

# exploring the variables in search of inconsistencies or typos

# species
table(sort(body.variance.full$species))


# trait
table(body.variance.full$trait) #in reality the all mean more or less the same, (flattened vs unflattened is the main difference), and we won't make a difference for the analyses.


# mean.male and mean.female
summary(body.variance.full$mean.male)
summary(body.variance.full$mean.female)

# checking for gross outliers
plot(body.variance.full$mean.male,body.variance.full$mean.female)
abline(1,1)

# only bustards are a bit of an outlier, but that makes sense
# body.variance.full[body.variance.full$mean.male>500,]


# n.male and n.female
summary(body.variance.full$n.male)
summary(body.variance.full$n.female)

# "extreme" values for Euplectes orix, from Craig et al. 2010; and also for Ficedula albicollis: double-checking the original references indicates that they are correct
# body.variance.full[body.variance.full$n.male>3000,]
# body.variance.full[body.variance.full$n.female>3000,]


# SD.male and SD.female
summary(body.variance.full$SD.male)
summary(body.variance.full$SD.female)
# checking SD == 0, it is suspicious
body.variance.full[body.variance.full$SD.male==0,] #sample size is very small, so by chance, the three individuals could have had indeed the same exact measurement. Leave it for now.


# do we have evidence for Taylor's law?
# males
plot(log(body.variance.full$mean.male),log(body.variance.full$SD.male))

# females
plot(log(body.variance.full$mean.female),log(body.variance.full$SD.female))


# reference
table(body.variance.full$reference)


# sorting database alphabetically
body.variance.full <- body.variance.full[order(body.variance.full$species),]

# creating an identifier for those male and female estimates that were at least reported together,
# since we are going to use female and male data as independent data points from now on
body.variance.full$pairID <- 1:nrow(body.variance.full)


# excluding estimates based on a sample size smaller than 2
# splitting the database so that males and females can be analyzed together
body.variance.full.final.males <- as.data.frame(body.variance.full 
                                                %>% select(-mean.female, -n.female, -SD.female) 
                                                %>% filter(n.male > 1 & !is.na(SD.male)) #excluding estimates based on less than 2 data points for the time being
                                                %>% rename(mean=mean.male,n=n.male,SD=SD.male)
                                                %>% mutate(sex="male"))

body.variance.full.final.females <- as.data.frame(body.variance.full 
                                                  %>% select(-mean.male, -n.male, -SD.male)
                                                  %>% filter(n.female > 1 & !is.na(SD.female)) #excluding estimates based on less than 2 data points for the time being
                                                  %>% rename(mean=mean.female,n=n.female,SD=SD.female)
                                                  %>% mutate(sex="female"))

# building the combined male-female database for the analyses
body.variance.full.sex.sep <- rbind(body.variance.full.final.males,body.variance.full.final.females)
body.variance.full.sex.sep$sex <- as.factor(body.variance.full.sex.sep$sex)

summary(body.variance.full.sex.sep)


######################################################################################################
# Phylogeny construction and species list fixing
######################################################################################################

# # obtaining dataframe listing the Open Tree identifiers potentially matching our list of species (be aware that this will take a few minutes. You can avoid waiting by loading the data below)
# study.species <- as.character(unique(body.variance.full.sex.sep$species))
# 
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
# Based on the checks above, we had to implement some changes to the species names:
# The remaining changes are typos that needed correction, update in species names or bringing subspecies up to species level so that everything is comparable
body.variance.full.sex.sep$species.updated <- recode(body.variance.full.sex.sep$species,
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
                                                     #"Porphyri porphyrio" = "Porphyrio porphyrio", #typo fixed
                                                     
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
                                                     "Cyanopica cyanus" = "Cyanopica cooki", #updating name based on population of origing
                                                     
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
                                                     #"Anas cyanoptera septentrionalium" = "Anas cyanoptera", # not brought up to species because this species is no longer monophyletic, which means that only by keeping it as subspecies we can get it's phylogentiec position
                                                     "Calonectris diomedea diomedea" = "Calonectris diomedea", # although it could probably stay as subspecies and nothing would change
                                                     "Emberiza schoeniclus lusitanica" = "Emberiza schoeniclus", # although it could probably stay as subspecies and nothing would change
                                                     "Emberiza schoeniclus schoeniclus" = "Emberiza schoeniclus", # although it could probably stay as subspecies and nothing would change
                                                     "Emberiza schoeniclus witherbyi" = "Emberiza schoeniclus", # although it could probably stay as subspecies and nothing would change
                                                     "Haematopus ostralegus longipes" = "Haematopus ostralegus", # although it could probably stay as subspecies and nothing would change
                                                     "Haematopus ostralegus ostralegus" = "Haematopus ostralegus", # although it could probably stay as subspecies and nothing would change
                                                     "Hydrobates pelagicus pelagicus" = "Hydrobates pelagicus", #brought down to species level because there are other instances with Hydrobates pelagicus only, just to make it comparable
                                                     "Malurus leucopterus edouardi" = "Malurus leucopterus", # brought down to species because subsp ranges largely overlap and also because this was the only species with multiple subspecies left, meaning that brining it down to the species level makes it more comparable to the rest of the dataset
                                                     "Malurus leucopterus leuconotus" = "Malurus leucopterus", # brought down to species because subsp ranges largely overlap and also because this was the only species with multiple subspecies left, meaning that brining it down to the species level makes it more comparable to the rest of the dataset
                                                     "Malurus leucopterus leucopterus" = "Malurus leucopterus", # brought down to species because subsp ranges largely overlap and also because this was the only species with multiple subspecies left, meaning that brining it down to the species level makes it more comparable to the rest of the dataset
                                                     "Zonotrichia leucophrys oriantha" = "Zonotrichia leucophrys", # although it could probably stay as subspecies and nothing would change
                                                     "Gallinula galeata sandvicensis" = "Gallinula galeata",
                                                     .default = levels(body.variance.full.sex.sep$species))

body.variance.full.sex.sep$species.updated <- factor(body.variance.full.sex.sep$species.updated)

# # repeating the process after fixing the species names
# # obtaining dataframe listing the Open Tree identifiers potentially matching our list of species (be aware that this will take a few minutes. You can avoid waiting by loading the data below)
# study.species <- as.character(unique(body.variance.full.sex.sep$species.updated))
# 
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
# # saving the taxonomic data created on the 7th of December 2021 to speed the process in the future
# save(taxa,file = "data/taxa_Open_Tree_of_Life_20211207.RData")

# loading the taxonomic data created on the 31st of August 2021
load("data/taxa_Open_Tree_of_Life_20211207.RData") #taxa


# everything seems in order now

# # retrieving phylogenetic relationships among taxa in the form of a trimmed sub-tree
# tree <- tol_induced_subtree(ott_ids = taxa[["ott_id"]], label_format = "name")
# 
# # we need to check for the existence of polytomies
# is.binary(tree) # there is a polytomy for Anous stolidus and Chroicocephalus ridibundus (see https://doi.org/10.1111/j.1469-7998.2009.00622.x)
# 
# # to take care of this polytomy, we are going to use a randomization approach
# set.seed(23)
# tree_random <- multi2di(tree,random=TRUE)
# is.binary(tree_random)
# 
# # to confirm that our tree covers all the species we wanted it to include, and make sure that the species names in our database match those in the tree, we use the following code
# 
# tree_random$tip.label <- gsub("_"," ", tree_random$tip.label)
# intersect(as.character(tree_random$tip.label), as.character(body.variance.full.sex.sep$species.updated))
# setdiff(as.character(body.variance.full.sex.sep$species.updated), as.character(tree_random$tip.label)) #listed in our database but not in the tree
# setdiff(as.character(tree_random$tip.label),as.character(body.variance.full.sex.sep$species.updated)) # listed in the tree but not in our database
# 
# 
# # we can now save the tree
# save(tree_random, file = "data/tree_random.Rdata")
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
# save(phylo_cor, file = "data/phylo_cor.Rdata")


# loading tree and matrix
load("data/tree_random.Rdata")
load("data/phylo_cor.Rdata")


# # we can then plot the tree
# plot(tree_random, type = "fan", cex=0.8, label.offset =.1, no.margin = TRUE) #check: https://www.rdocumentation.org/packages/ape/versions/5.3/topics/plot.phylo
# 
# plot(tree_random, cex=0.25, label.offset =.1, no.margin = TRUE) #check: https://www.rdocumentation.org/packages/ape/versions/5.3/topics/plot.phylo


######################################################################################################
# Adding migratory and feeding type information
######################################################################################################

# first creating a database to be exported and filled in with data
species.database <- unique(body.variance.full.sex.sep[,c("species.updated","location")])

# adding variables needed to be filled up
species.database$migratory <- ""
species.database$migratory.comments <- ""
species.database$feeding.type.aerial <- ""
species.database$feeding.type.comments <- ""

species.database <- unique(species.database[,c("species.updated","migratory","migratory.comments",
                                               "feeding.type.aerial","feeding.type.comments","location")])

species.database <- species.database[order(species.database$species.updated),]

# saving database so that AST and mostly KR can fill this information up for the analyses
# write.csv(unique(species.database[,c("species.updated","migratory","migratory.comments",
#                                      "feeding.type.aerial","feeding.type.comments")]),
#           "data/final/final_bird_species_list.csv",row.names=FALSE)


######################################################################################################
# Then, adding the species-specific information about migration type and feeding type

# importing data
species.information <- read.csv("data/final/final_bird_species_list_migration_and_feeding.csv",header=T,sep=",")

# reformatting some variables types to factors
cols.factors.2 <- c("species.updated","migratory","feeding.type.aerial")
species.information[cols.factors.2] <- lapply(species.information[cols.factors.2], factor)


# merging databases to create the final database for analyses
body.variance.full.sex.sep.species.info <- merge(body.variance.full.sex.sep,species.information,
                                                 all.x = T, by = "species.updated")


# last fix to the dataset, removing an estimate based on three individuals which has an SD==0
body.variance.full.sex.sep.species.info <- as.data.frame(body.variance.full.sex.sep.species.info %>% filter(SD>0))


######################################################################################################
# Calculating lnCV
######################################################################################################

# calculating lnCV
body.variance.full.sex.sep.species.info <- as.data.frame(escalc(measure = "CVLN", mi=mean, sdi=SD, ni=n,
                                                                data=body.variance.full.sex.sep.species.info,
                                                                add.measure=FALSE,append=TRUE))


# generating a new species variable to account for repeated values for some species + phylogeny
body.variance.full.sex.sep.species.info$species.updated.rep <- body.variance.full.sex.sep.species.info$species.updated

# generating a unit level effect to model residual/within-study variance
body.variance.full.sex.sep.species.info$effectsizeID <- 1:nrow(body.variance.full.sex.sep.species.info)


# # quick funnel plots
# # precision
# plot(sqrt(1/body.variance.full.sex.sep.species.info$vi),exp(body.variance.full.sex.sep.species.info$yi))
# abline(exp(mean(body.variance.full.sex.sep.species.info$yi,na.rm=T)),0)
# abline(exp(median(body.variance.full.sex.sep.species.info$yi,na.rm=T)),0,col="red")
# 
# #se
# plot(sqrt(body.variance.full.sex.sep.species.info$vi),exp(body.variance.full.sex.sep.species.info$yi))
# abline(exp(mean(body.variance.full.sex.sep.species.info$yi,na.rm=T)),0)
# abline(exp(median(body.variance.full.sex.sep.species.info$yi,na.rm=T)),0,col="red")
# 
# #se
# plot(body.variance.full.sex.sep.species.info$vi,exp(body.variance.full.sex.sep.species.info$yi))
# abline(exp(mean(body.variance.full.sex.sep.species.info$yi,na.rm=T)),0)
# abline(exp(median(body.variance.full.sex.sep.species.info$yi,na.rm=T)),0,col="red")
# 
# #n original scale CV
# plot(body.variance.full.sex.sep.species.info$n,(body.variance.full.sex.sep.species.info$yi))
# abline((mean(body.variance.full.sex.sep.species.info$yi,na.rm=T)),0)
# abline((median(body.variance.full.sex.sep.species.info$yi,na.rm=T)),0,col="red")
# 
# #n
# plot(body.variance.full.sex.sep.species.info$n,exp(body.variance.full.sex.sep.species.info$yi))
# abline(exp(mean(body.variance.full.sex.sep.species.info$yi,na.rm=T)),0)
# abline(exp(median(body.variance.full.sex.sep.species.info$yi,na.rm=T)),0,col="red")
# 
# 
# plot(body.variance.full.sex.sep.species.info$n,exp(body.variance.full.sex.sep.species.info$yi),xlim=c(0,100),ylim=c(0,0.16))

# there are a couple of suspiciously high values, but checking the original papers did not clarify whether they are mistakes (SE vs SD), typos or just real, therefore, we are leacing them in


