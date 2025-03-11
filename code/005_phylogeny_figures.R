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

# This script is to generate Figure 1 in the main manuscript for the following
# study:

# Klaus Reinhold, Alfredo Sánchez-Tójar. 2025. Wing length canalisation and 
# behaviour across birds: a phylogenetic meta-analysis of variance

################################################################################
# Packages needed
################################################################################

# install.packages("pacman")
pacman::p_load(metafor,ggplot2,ggcorrplot,
               rgbif,tidyverse,
               ape,rotl,ggtree,rsvg)

# cleaning environment
rm(list=ls())

################################################################################
# Importing the data 
################################################################################

# importing final dataset
body.CV.final <- read.csv("data/final/04_final_full_and_clean_bird_wing_size_dataset.csv",
                          header=T,sep=",")

# quick view
names(body.CV.final)
summary(body.CV.final)


# converting some variables types to factors
cols.factors <- c("population_ID","species.updated","species",
                  "study_ID","reference_link","subset",
                  "trait","unit","measurement_notes",
                  "pair_ID","sex","migratory","feeding.type",
                  "migratory.pop.level","feeding.type.pop.level",
                  "migratory.inc.vagrants.pop.level","migratory.inc.vagrants")

body.CV.final[cols.factors] <- lapply(body.CV.final[cols.factors], 
                                      factor)

summary(body.CV.final)

# adding the missing random effect variables needed for the analyses

# generating a new species variable to account for repeated values for some 
# species + phylogeny (sensu Cinar et al. 2020, MEE)
body.CV.final$species.updated.rep <- body.CV.final$species.updated

# generating a unit level effect to model residual/within-study variance
body.CV.final$effectsize_ID <- 1:nrow(body.CV.final)


################################################################################
# Calculating our effect size of interest: lnCV
################################################################################

# calculating lnCV. Note that escalc() includes a small-size bias correction
body.CV.final <- escalc(measure = "CVLN",
                        mi=mean,
                        sdi=SD,
                        ni=n,
                        data=body.CV.final,
                        add.measure=FALSE,
                        append=TRUE)

# back-transforming lnCV to CV (coefficient of variation) to aid visualization
body.CV.final$CV <- exp(body.CV.final$yi)*100

summary(body.CV.final)

################################################################################
# Loading phylogenetic information at the species level

# taxonomic information: taxa
load("data/phylogeny/taxa_Open_Tree_of_Life_20240909.RData")

# tree: tree_random
load("data/phylogeny/tree_random_20240909.Rdata")

# phylogenetic matrix: phylo_cor
load("data/phylogeny/phylo_cor_20240909.Rdata")


################################################################################
# Adding extra taxonomic information: genus, family, order
################################################################################

# adding extra taxonomic information
name_list <- unique(as.character(body.CV.final$species.updated))

# code kindly provided by Matt Grainger (Thank you!)
# taxonomic.info.full <- name_backbone_checklist(name_list)# |> view()
# 
# # saving the taxonomic data to speed the process in the future and allow full
# # computational reproducibility
# save(taxonomic.info.full,file = "data/phylogeny/taxonomic_extra_information_20241109.RData")
# loading the most updated taxonomic data created
load("data/phylogeny/taxonomic_extra_information_20241109.RData") #taxonomic.info.full

# extracting our information of interest
taxonomic.info.full <- as.data.frame(taxonomic.info.full)[,c("verbatim_name",
                                                             "species",
                                                             "genus",
                                                             "order",
                                                             "family")]
# # let's check it out
# taxonomic.info.full

# making sure there is full overlap between the two datasets and that we are
# using the correct variable to merge them
setdiff(name_list,
        as.character(taxonomic.info.full$verbatim_name))

setdiff(as.character(taxonomic.info.full$verbatim_name),
        name_list)


# manually completing the missing information for greenfinches
taxonomic.info.full[taxonomic.info.full$verbatim_name=="Chloris chloris",] <- 
  c("Chloris chloris","Chloris chloris","Chloris","Passeriformes","Fringillidae")

# merging the datasets
taxonomic.info.full <- taxonomic.info.full[order(taxonomic.info.full$order),]
taxonomic.info.full

# adding taxonomic information to our database
body.CV.final <- merge(body.CV.final,
                       taxonomic.info.full,
                       by.x = "species.updated",
                       by.y = "verbatim_name")

summary(body.CV.final)

# removing species.y
body.CV.final <- dplyr::select(body.CV.final, -species.y)

################################################################################
# how many species, genera, families and orders do we have?
length(unique(body.CV.final$species.updated))
length(unique(body.CV.final$genus))
length(unique(body.CV.final$family))
length(unique(body.CV.final$order))

# what's the number of species per order?
as.data.frame(unique(body.CV.final[,c("species.updated","family","order")]) %>%
                group_by(order) %>%
                summarise(N_obs = n()))  # calculate sample size

# what's the number of species per family?
as.data.frame(unique(body.CV.final[,c("species.updated","family","order")]) %>%
                group_by(family) %>%
                summarise(N_obs = n()))  # calculate sample size


################################################################################
# Full phylogenetic tree
################################################################################

# Resources use to generate this phylogeny figure:
# https://yulab-smu.top/treedata-book/index.html
# https://stackoverflow.com/questions/44783855/ggtree-plotting-area-not-big-enough
# https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html

# first, generating a dataset with mean CV-values per species for visualizing 
# purposes using clustering
agg.es.species <- as.data.frame(escalc(measure = "CVLN", 
                                       mi = mean, 
                                       sdi = SD, 
                                       ni = n,
                                       data = body.CV.final) %>%
                                  aggregate(cluster = species.updated, # using species.updated as clustering variable
                                            struct = "CS", # specifying varcovar compound symmetric structure of the sampling errors
                                            rho = 0.5, # with an assumed correlation of 0.5
                                            addk = TRUE))

agg.es.species$CV <- exp(agg.es.species$yi)*100

# first basic tree to build upon
full.tree <- ggtree(tree_random, 
                    layout='circular',
                    size = 1.05)# + 
#ggtitle("Phylogenetic tree (species-level)")

full.tree


#fulltree.1 <- full.tree %<+% body.CV.final[,c("species.updated","order","CV")]
fulltree.1 <- full.tree %<+% agg.es.species[,c("species.updated","order","CV")]
fulltree.1

fulltree.2 <- fulltree.1 + 
  geom_tippoint(aes(color = order, size = CV)) + 
  scale_size(range=c(0.5,3)) +
  geom_tiplab(size = 2.2, offset = .9, aes(angle = angle)) +
  guides(colour = guide_legend(title = "Order"),
         size = guide_legend(title = "Wing length CV (%)")) +
  theme(legend.position = c(1, 0.52),
        #legend.position = "right",
        #legend.justification = c(1, 0.5),
        plot.margin = unit(c(t=-2,r=-4,b=-2,l=-8), "cm"),
        legend.key.size = unit(0.55, units = "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10, face = "bold"),
        title = element_text(size = 15, face = "bold"))# + 
#scale_size_continuous(range = c(3, 10))
fulltree.2

# checking nodes
tree_random$tip.label

# generating a list of node positions where to put the silhouettes
node.list <- c(11,22,28,43,55,68,
               75,87,97,102,106,122,
               132,142,144,152,159,
               161,168)
# here is the unique ids (from the phylopic website link) for each node
image.list <- c("ec11f62e-ed3e-49bb-bcff-4e47b281d378", # turdus merula node: 11
                "19f3f55c-f942-464b-a61f-27794b9000f7", # geospiza fortis node: 20
                "fb0ae249-66cc-4fce-b668-22ffabcb6bdc", # emberiza node: 28
                "3a4cdd72-e553-40ad-838e-3b23037b2010", # passer domesticus node: 43
                "49d94179-c795-4a53-938b-9c4d2d2ac692", # aegithalos caudatus node: 55
                "784ed1b5-0217-421c-8397-289c7cbb472a", # aludinae node: 68
                "715f6434-de54-475b-9653-4b4db6864bdc", # cyanopica cooki node: 73
                "ddd5783c-ded5-48f2-a07d-cc37c83b227b", # dendrocopos node: 90
                "63f337b6-8a58-4939-ba8b-f1bd3b8c1037", # milvus migrans node: 93
                "ef6d1b7a-15fe-46cc-aa1e-2598d767895c", # calidris node: 102
                "b99b9cef-88bf-455f-8190-675aa11c1fb2", # calidris pugnax node: 106
                "18ff6244-3ec9-4750-a29a-aaad7e2f14ad", # fratercula arctica node: 122
                "deba1d91-daa8-40a6-8d48-7a9f295bc662", # podiceps node: 132
                "56554944-1a1c-4281-a6e5-49ce2de9b2f4", # pterodroma node: 142
                "9171cd2b-3afc-46f4-9ee1-c6515de0378c", # phaethon rubricauda node: 145
                "cd5ca1a2-b0a9-4ce9-b163-5c004e327d9d", # rallus aquaticus node: 152
                "737a0176-deee-4a0f-a721-11534ff2b0c5", # apus pallidus node: 159
                "42f85a2b-7517-439a-8dc3-b745a35c035d", # lagopus lagopus node: 161
                "ad8b9af1-d189-4fb7-b555-fcd7f27afccb") # branta leucopsis node: 168

# generating a dataset for adding the images
dt <- data.frame(node = node.list,
                 name = rep("sp",length(node.list)),
                 image = image.list) 
fulltree.3 <- fulltree.2 + geom_cladelab(data = dt, 
                                         mapping = aes(node = node, 
                                                       label = name, 
                                                       image = image), 
                                         geom = "phylopic", 
                                         imagecolor = "black", 
                                         offset=20, offset.text=-3)

# saving the figure
png(filename = 'figures/Figure_1_Panel_1_Phylogenetic_tree_Species.png', 
    width = 27, height = 21, units = 'cm', 
    res = 600)

fulltree.3

dev.off()


################################################################################
# Visualizing the phylogenetic correlation matrix
################################################################################

# from https://yefeng0920.github.io/MA_Map_Bib/#phylogenetic-trees

# saving the figure
png(filename = 'figures/Figure_SX_Phylogenetic_varcovar_correlation_matrix.png', 
    width = 30, height = 30, units = 'cm', 
    res = 800)

# visual
ggcorrplot::ggcorrplot(phylo_cor, sig.level = 0.05, lab_size = 1,
                       p.mat = NULL,insig = c("pch", "blank"),
                       pch = 1, pch.col = "black", pch.cex = 1, tl.cex = 1.5) +
  theme(axis.text.x = element_text(size = 4, margin = margin(-2, 0, 0, 0)),
        axis.text.y = element_text(size = 4, margin = margin(0, -2, 0, 0)),
        panel.grid.minor = element_line(size = 3)) +
  geom_tile(fill = "white") +
  geom_tile(height = 0.8, width = 0.8) +
  scale_fill_gradient2(low = "#E69F00",mid = "white", high = "#56B4E9",
                       midpoint = 0.5, breaks = c(0, 1),
                       limit = c(0,1)) + labs(fill = "Correlation")

dev.off()

################################################################################
# Building a tree at the Order level for visualizing order mean CV values
################################################################################

# # extracting taxonomic information
# taxa.order <- tnrs_match_names(names = unique(body.CV.final$order))
# taxa.order # some orders seem already problematic: barrent, merged
# 
# # check whether these orders occur in the synthetic tree
# ott_in_tree <- ott_id(taxa.order)[is_in_tree(ott_id(taxa.order))]
# length(ott_id(taxa.order)) - length(is.na(ott_in_tree)) 
# # There are 6 orders that don't appear in the synthetic tree
# # (doing this at the family level would end up with 39 missing families)
# 
# ott_not_in_tree <- ott_id(taxa.order)[!is_in_tree(ott_id(taxa.order))]
# is_in_tree(ott_not_in_tree)
# # > is_in_tree(ott_not_in_tree)
# # Charadriiformes Gruiformes Coraciiformes Otidiformes Phaethontiformes Suliformes 
# # FALSE           FALSE      FALSE         FALSE       FALSE            FALSE 
# ott_not_in_tree

# Since there are 6 orders that do not appear in the synthetic tree but we 
# still want to show a tree with the 16 orders so that we can visualize 
# effect sizes and orders visually, we are going to trick the tree for that
# purpose by manually changing the names of orders that were not found to their
# corresponding family names

body.CV.final$order.fake <- ifelse(body.CV.final$order=="Charadriiformes",
                                   #"Charadriidae",#charadriidae, laridae are also not in the tree
                                   "Scolopacidae",
                                   body.CV.final$order)

body.CV.final$order.fake <- ifelse(body.CV.final$order.fake=="Gruiformes",
                                   "Gruidae",
                                   body.CV.final$order.fake)

body.CV.final$order.fake <- ifelse(body.CV.final$order.fake=="Coraciiformes",
                                   "Coraciidae",
                                   body.CV.final$order.fake)

body.CV.final$order.fake <- ifelse(body.CV.final$order.fake=="Otidiformes",
                                   "Otididae",
                                   body.CV.final$order.fake)

body.CV.final$order.fake <- ifelse(body.CV.final$order.fake=="Phaethontiformes",
                                   "Phaethontidae",
                                   body.CV.final$order.fake)

body.CV.final$order.fake <- ifelse(body.CV.final$order.fake=="Suliformes",
                                   "Sulidae",
                                   body.CV.final$order.fake)

unique(body.CV.final$order.fake)


################################################################################
# building family/order tree
# Extra information here: https://lunasare.github.io/ssb2020_workshop/01-find-taxa/index.html
# taxa.order <- tnrs_match_names(names = unique(body.CV.final$order.fake))
# taxa.order
# 
# # check whether occur in the synthetic tree: all in there now
# ott_in_tree <- ott_id(taxa.order)[is_in_tree(ott_id(taxa.order))]
# length(ott_id(taxa.order)) - length(is.na(ott_in_tree))
# 
# 
# # saving the taxonomic data to speed the process in the future and allow full 
# # computational reproducibility
# save(taxa.order,file = "data/phylogeny/taxa_Open_Tree_of_Life_Order_20241109.RData")
# # loading the most updated taxonomic data created
load("data/phylogeny/taxa_Open_Tree_of_Life_Order_20241109.RData") #taxa.order


# # retrieving phylogenetic relationships among taxa in the form of a trimmed
# # sub-tree
# tree.order <- tol_induced_subtree(ott_ids = taxa.order[["ott_id"]])
# 
# # To make sure the tree combining both orders and families makes sense, I went
# # to https://tree.opentreeoflife.org/opentree/opentree15.1@ott81461/Aves
# # and make sure things align, which they did.
# 
# # we need to check for the existence of polytomies
# is.binary(tree.order)
# 
# # making sure both tree and dataset agree
# # here are the species included in the tree
# sort(tree.order$tip.label)
# 
# # removing the underscore "_" and ott id from the tree tip.label
# tree.order$tip.label <- str_split(tree.order$tip.label, "_", simplify = TRUE)[,1]
# sort(tree.order$tip.label)
# 
# # here are the species included in our dataset
# sort(unique(as.character(body.CV.final$order.fake)))
# 
# # comparing the tree to our dataset
# intersect(as.character(tree.order$tip.label),
#           as.character(body.CV.final$order.fake))
# setdiff(as.character(body.CV.final$order.fake),
#         as.character(tree.order$tip.label)) # listed in our database but not in the tree
# setdiff(as.character(tree.order$tip.label),
#         as.character(body.CV.final$order.fake)) # listed in the tree but not in our database
# 
# # compute branch lengths of tree
# phylo_branch_order <- compute.brlen(tree.order, method = "Grafen", power = 1)
# 
# # we can now save the tree
# save(phylo_branch_order, 
#      file = "data/phylogeny/tree_branch_Order_20241109.RData")
# loading tree
load("data/phylogeny/tree_branch_Order_20241109.RData") #phylo_branch_order

# check tree is ultrametric
is.ultrametric(phylo_branch_order) # TRUE

# matrix to be included in the models
phylo_cor_order <- vcv(phylo_branch_order, cor = T)

# visually exploring the vcv for this tree at the order level
ggcorrplot::ggcorrplot(phylo_cor_order, sig.level = 0.05, lab_size = 1,
                       p.mat = NULL,insig = c("pch", "blank"),
                       pch = 1, pch.col = "black", pch.cex = 1, tl.cex = 1.5) +
  theme(axis.text.x = element_text(size = 4, margin = margin(-2, 0, 0, 0)),
        axis.text.y = element_text(size = 4, margin = margin(0, -2, 0, 0)),
        panel.grid.minor = element_line(size = 3)) +
  geom_tile(fill = "white") +
  geom_tile(height = 0.8, width = 0.8) +
  scale_fill_gradient2(low = "#E69F00",mid = "white", high = "#56B4E9",
                       midpoint = 0.5, breaks = c(0, 1),
                       limit = c(0,1)) + labs(fill = "Correlation")

# check tree labels
phylo_branch_order$tip.label

# changing the family names to their original corresponding order names
phylo_branch_order$tip.label <- ifelse(phylo_branch_order$tip.label=="Scolopacidae",
                                       "Charadriiformes",
                                       phylo_branch_order$tip.label)

phylo_branch_order$tip.label <- ifelse(phylo_branch_order$tip.label=="Gruidae",
                                       "Gruiformes",
                                       phylo_branch_order$tip.label)

phylo_branch_order$tip.label <- ifelse(phylo_branch_order$tip.label=="Coraciidae",
                                       "Coraciiformes",
                                       phylo_branch_order$tip.label)

phylo_branch_order$tip.label <- ifelse(phylo_branch_order$tip.label=="Otididae",
                                       "Otidiformes",
                                       phylo_branch_order$tip.label)

phylo_branch_order$tip.label <- ifelse(phylo_branch_order$tip.label=="Phaethontidae",
                                       "Phaethontiformes",
                                       phylo_branch_order$tip.label)

phylo_branch_order$tip.label <- ifelse(phylo_branch_order$tip.label=="Sulidae",
                                       "Suliformes",
                                       phylo_branch_order$tip.label)

# check tree labels
phylo_branch_order$tip.label


################################################################################
# Plotting taxonomic effect sizes
################################################################################

# following code provided at: https://yefeng0920.github.io/MA_Map_Bib/#phylogenetic-trees

# calculating CV values for each order using clustering
agg.es <- escalc(measure = "CVLN", 
                 mi = mean, 
                 sdi = SD, 
                 ni = n,
                 data = body.CV.final) %>%
  aggregate(cluster = order, # using order as clustering variable
            struct = "CS", # specifying varcovar compound symmetric structure of the sampling errors
            rho = 0.5, # with an assumed correlation of 0.5
            addk = TRUE)

# subsetting and renaming the dataframes
agg.es <- agg.es[c("order", "yi", "vi","ki")]
names(agg.es) <- c("Order", "Mean", "SE","N_obs")

# calculating 95% CI for each order, keep in mind that this calculation is for
# exploration purposes since no random effects are accounted for there
agg.es <- agg.es %>%
  mutate(Lower_bound = Mean - sqrt(SE) * qnorm(0.975), 
         Upper_bound = Mean + sqrt(SE) *
           qnorm(0.975)) %>%
  arrange(Order)

agg.es

# calculating the number of observations per order
N_obs_k <- body.CV.final %>%
  group_by(order) %>%
  summarise(N_obs = n())  # calculate sample size

# calculating the number of species per order
N_obs <- body.CV.final %>%
  group_by(order) %>%
  summarise(N_species = n_distinct(species.updated))

# renaming variable
as.data.frame(N_obs[order(N_obs$N_species),])
names(N_obs)[1] <- "Order"

# adding number of species per order to ourdataset for the plot
fe.spp.es <- left_join(agg.es, N_obs, by = "Order")

# tip.label <- data.frame(Order = phylo_branch_order$tip.label)  # extract tip label
# fe.spp.es2 <- left_join(tip.label, fe.spp.es, by = "Order")

# calculaing z values and p-values for those estimates, which in our case is 
# kind of pointless because they are all very different from zero anyway, and 
# we do no need to plot those z and p-values
fe.spp.es2 <- fe.spp.es %>%
  mutate(z = Mean/SE, p = pnorm(abs(z), lower.tail = F) * 2)

# backtransforming lnCV to CV% to aid biological interpretation
fe.spp.es2$Mean.perc <- exp(fe.spp.es2$Mean)*100
fe.spp.es2$Lower_bound.perc <- exp(fe.spp.es2$Lower_bound)*100
fe.spp.es2$Upper_bound.perc <- exp(fe.spp.es2$Upper_bound)*100

# here is our final dataset for the order-level plots
fe.spp.es2

# make the tree by adding plot elements sequentially:
tree.p1 <- ggtree(phylo_branch_order, layout = "rectangular", cex = 0.4) 
tree.p1

tree.p1.2 <- tree.p1 %<+% as.data.frame(fe.spp.es2[,c("Order","Order","Order",
                                                      "N_obs","N_obs")]) #trick
tree.p1.2


tree.p2 <- tree.p1.2 +
  geom_tiplab(aes(color = Order.2), 
              #color="black",
              size = 5, fontface = "italic",
              align = T, offset = 0.05) +
  geom_tippoint(aes(size = N_obs,color = Order.1)) +
  #guides(color = "none") +
  theme(legend.position="none") +
  guides() +
  xlim_expand(xlim = c(0, 1.5), panel = "Tree") 

tree.p2

tree.p3 <- facet_plot(tree.p2, 
                      panel = "Number of species", 
                      data = fe.spp.es2, 
                      geom = ggstance::geom_barh,
                      #mapping = aes(x = N_obs),
                      mapping = aes(x = N_species,  # showing the number of species
                                    fill = Order.1),
                      #fill = "lightblue", color = "blue4",
                      color = "blue4",
                      #size = "none",
                      alpha = 0.5, stat = "identity") +
  guides(size = "none")

tree.p3

tree.p4 <- tree.p3 + 
  geom_facet(panel = "Wing length CV (%)",
             data = fe.spp.es2,
             geom = ggstance::geom_pointrangeh,
             mapping = aes(x = Mean.perc,
                           xmin = Lower_bound.perc,
                           xmax = Upper_bound.perc,
                           color = Order.1,
                           size = N_obs.1/4))+
  scale_size(range=c(0.75,3.25)) +
  theme_tree2() + 
  theme(strip.background = element_rect(fill = "lightblue"),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold"
        )) + 
  guides(fill = "none",
         color = "none") +
  xlim_expand(xlim = c(1.5, 6), panel = "Wing length CV (%)")


tree.p4

tree.p5 <- facet_widths(tree.p4, 
                        c(`Tree` = 0.5, 
                          `Number of species` = 0.2, 
                          `Wing length CV (%)` = 0.3))  #adjust widths of the plot facets for display


# saving the figure
png(filename = 'figures/Figure_1_Panel_2_Phylogenetic_tree_Order.png', 
    width = 27, height = 10, units = 'cm', 
    res = 600)

tree.p5

dev.off()



# plot_list(fulltree.3,
#           tree.p5,
#           align = "v", 
#           nrow = 2, 
#           rel_heights = c(3/4, 1/4),
#           ncol = 1,
#           tag_levels = c("A","B"))


# library(ggpubr)
# figure <- ggarrange(fulltree.3, 
#                     tree.p5,
#                     labels = c("A", "B"),
#                     ncol = 1, nrow = 2,
#                     heights=c(3,1))
# figure
# 
# # ggexport(bxp, dp, lp, bxp, filename = "test.pdf",
# #          nrow = 2, ncol = 1)