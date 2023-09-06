# another epidemic overlap method
  # instead of taking the minimum of each prevalance 
  # (or inf host density) across hosts, average the values
  # if both hosts were infected that day

rm(list=ls())

setwd("~/Dropbox (University of Michigan)/Manuscripts/McIntireetal_Networks/Data&Code/daphnia-parasite-networks/2023/networks")

# NOTE: requires all of the functions in "centrality-lake-functions.R"
source("centrality-lake-functions.R")

# packages
require(tidyverse)
require(zoo)
require(igraph)
require(bipartite)
require(dunn.test)
require(conover.test)
require(cowplot)


clean_data <- read.csv("Clean-Data-2014-2016.csv", header = TRUE)

# STEP 1: prep data, epidemic overlap of prevalence
# lakeyear_auc_shared_prev <- read.csv("lakeyear_auc_shared_prev.csv", header = TRUE)
lakeyear_auc_density <- read.csv("lakeyear_auc_density.csv", header = TRUE)
#lakeyear_auc_shared_prevdens <- read.csv("lakeyear_auc_shared_prevdens.csv", header = TRUE)
lakeyear_auc_shared_total_density <- read.csv("lakeyear_auc_shared_total_density.csv", header = TRUE)


# General code ---------------------------------------------------------------------
clean_auc <- read.csv("Clean-Data-2014-2016.csv", header = TRUE)

  # rename columns; easier to "gather" this way
  colnames(clean_auc)[which(names(clean_auc) == "pasteuria.prev")] <- "pasteuria"
  colnames(clean_auc)[which(names(clean_auc) == "metsch.prev")] <- "metsch"
  colnames(clean_auc)[which(names(clean_auc) == "scarlet.prev")] <- "scarlet"
  colnames(clean_auc)[which(names(clean_auc) == "brood.prev")] <- "brood"
  colnames(clean_auc)[which(names(clean_auc) == "larssonia.prev")] <- "larssonia"
  colnames(clean_auc)[which(names(clean_auc) == "spider.prev")] <- "spider"
  colnames(clean_auc)[which(names(clean_auc) == "gurleya.prev")] <- "gurleya"
  
  # split data depending on total number of hosts counted in parasite sample
  total_less <- clean_auc %>%
    filter(Total < 20)
  
  # if less than 20 total, set prevalence to 0
  total_less$pasteuria <- 0
  total_less$metsch <- 0
  total_less$scarlet <- 0
  total_less$brood <- 0
  total_less$larssonia <- 0
  total_less$spider <- 0
  total_less$gurleya <- 0
  
  # don't do anything if total is greater/equal to 20
  total_greater <- clean_auc %>%
    filter(Total >= 20)
  
  # join all data back together
  loop_data <- full_join(total_less, total_greater)
  
  # gather into tall format; select appropriate columns
  loop_data <- loop_data %>%
    gather(Parasite.Species, Prevalence, pasteuria:gurleya) %>%
    select(Host.Species, Year, Lake, Parasite.Species, Julian.Day, Prevalence)
  
  # arrange data 
  loop_data <- loop_data %>%
    group_by(Host.Species, Year, Lake, Parasite.Species, Julian.Day) %>%
    arrange()
  
  # set year and julian day as correct class
  loop_data$Year <- as.factor(loop_data$Year)
  loop_data$Julian.Day <- as.integer(loop_data$Julian.Day)
  
  # order everything by julian day
  loop_data <- loop_data[order(loop_data$Julian.Day), ] #make sure julian days are in order
  
  # make sure it's in data frame format
  loop_data <- as.data.frame(loop_data)
  
  
  
  
  
# Shared AUC loop (arithmetic mean) -------------------------------------------------------------------------
# calculates the overlapping prevalence area (same parasite, different hosts)

YEAR<-unique(loop_data$Year)

OUT <- NULL

for (i in YEAR) {
  thisyear <- loop_data[loop_data$Year == i, ] # look at 1 year at a time
  LAKE <- unique(thisyear$Lake)
  
  for (j in LAKE) {
    thislake <- thisyear[thisyear$Lake == j, ] # within a year, look at 1 lake at a time
    PARA <- unique(thislake$Parasite.Species)
    
    for (p in PARA) {
      thisparasite <- thislake[thislake$Parasite.Species == p, ] # within a lake, look at 1 parasite at a time
      SPECIES <- unique(thisparasite$Host.Species)
      
      for (m in SPECIES) { #where m = dent, retro, cerio, etc
        thisspecies <- thisparasite[thisparasite$Host.Species == m, ]
        stuffone <- select(thisspecies, Julian.Day, Prevalence) # focal host data
        
        for (n in SPECIES) {
          otherspecies <- thisparasite[thisparasite$Host.Species == n, ]
          stufftwo <- select(otherspecies, Julian.Day, Prevalence) # other host data
          both <- full_join(stuffone, stufftwo, by = "Julian.Day") # use whatever merging function makes sense
          
          
          both$Prevalence.x[both$Prevalence.x == 0] <- NA  
          both$Prevalence.y[both$Prevalence.y == 0] <- NA  
          avg_epi <- rowMeans(both[, c("Prevalence.x", "Prevalence.y")], na.rm = FALSE)
          both$avg_epi <- avg_epi   
          both[is.na(both)] <- 0
          avgboth <- both
          # avgboth will have columns: Julian.Day, Prevalence.x, Prevalence.y, avg_epi
          TOTS <- NULL #make an empty matrix to include sequential prevalences. This will be used to calculate integrated areas.
          
          # now calculate the auc using the avgboth data
          for (q in 1:length(avgboth$Prevalence.x)) { #CHECK THAT THIS PART MAKES SENSE
            total <- avgboth[q, 4] # average prevalence prevalence
            jul <- avgboth[q, 1] # julian day
            tots <- c(q, jul, total) # new row in dataframe with q, julian day, and total
            TOTS <- rbind(TOTS, tots) # bind all rows together (different one for each lake/year/parasite)
            colnames(TOTS) <- c("q", "date", "total")
          }
          
          rownames(TOTS) <- NULL
          TOTS <- as.data.frame(TOTS)
          CURVE <- NULL #make an empty matrix for calulating each chunk of the integrated area by trapezoid rule.
          
          for (r in 1:(length(TOTS$total) - 1)) {
            curve <- 0.5*((TOTS[r, 3] + TOTS[r+1, 3])*(TOTS[r+1, 2] - TOTS[r, 2])) #trapezoid rule; this calculates area of each chunk
            CURVE <- c(CURVE, curve)
          } #write down each chunk
          
          area <- sum(CURVE) #take the sum
          output <- c(i, j, p, m, n, area)
          OUT <- rbind(OUT, output)
          
        }
      }
    }
  }
}




colnames(OUT) <- c("Year", "Lake", "Parasite.Species", "Focal.Host", "Other.Host", "AUC.prev")
# write.table(OUT, "auc_14to16_prev.txt", col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE) 

# tidy up the data frame
auc_shared <- as.data.frame(OUT)
row.names(auc_shared) <- c()
auc_shared$AUC.prev <- as.numeric(as.character(auc_shared$AUC.prev))

lakeyear_auc_shared_prev <- auc_shared


# Shared AUC loop (GEOMETRIC mean) -------------------------------------------------------------------------
# calculates the overlapping prevalence area (same parasite, different hosts)

YEAR<-unique(loop_data$Year)

OUT <- NULL

for (i in YEAR) {
  thisyear <- loop_data[loop_data$Year == i, ] # look at 1 year at a time
  LAKE <- unique(thisyear$Lake)
  
  for (j in LAKE) {
    thislake <- thisyear[thisyear$Lake == j, ] # within a year, look at 1 lake at a time
    PARA <- unique(thislake$Parasite.Species)
    
    for (p in PARA) {
      thisparasite <- thislake[thislake$Parasite.Species == p, ] # within a lake, look at 1 parasite at a time
      SPECIES <- unique(thisparasite$Host.Species)
      
      for (m in SPECIES) { #where m = dent, retro, cerio, etc
        thisspecies <- thisparasite[thisparasite$Host.Species == m, ]
        stuffone <- select(thisspecies, Julian.Day, Prevalence) # focal host data
        
        for (n in SPECIES) {
          otherspecies <- thisparasite[thisparasite$Host.Species == n, ]
          stufftwo <- select(otherspecies, Julian.Day, Prevalence) # other host data
          both <- full_join(stuffone, stufftwo, by = "Julian.Day") # use whatever merging function makes sense
          
          
          both$Prevalence.x[both$Prevalence.x == 0] <- NA  
          both$Prevalence.y[both$Prevalence.y == 0] <- NA  
          avg_epi <- apply(both[, c("Prevalence.x", "Prevalence.y")], 1, function(x) (prod(x[x!=0]))^(1/sum(x!=0)))
          both$avg_epi <- avg_epi   
          both[is.na(both)] <- 0
          avgboth <- both
          # avgboth will have columns: Julian.Day, Prevalence.x, Prevalence.y, avg_epi
          TOTS <- NULL #make an empty matrix to include sequential prevalences. This will be used to calculate integrated areas.
          
          # now calculate the auc using the avgboth data
          for (q in 1:length(avgboth$Prevalence.x)) { #CHECK THAT THIS PART MAKES SENSE
            total <- avgboth[q, 4] # average prevalence prevalence
            jul <- avgboth[q, 1] # julian day
            tots <- c(q, jul, total) # new row in dataframe with q, julian day, and total
            TOTS <- rbind(TOTS, tots) # bind all rows together (different one for each lake/year/parasite)
            colnames(TOTS) <- c("q", "date", "total")
          }
          
          rownames(TOTS) <- NULL
          TOTS <- as.data.frame(TOTS)
          CURVE <- NULL #make an empty matrix for calulating each chunk of the integrated area by trapezoid rule.
          
          for (r in 1:(length(TOTS$total) - 1)) {
            curve <- 0.5*((TOTS[r, 3] + TOTS[r+1, 3])*(TOTS[r+1, 2] - TOTS[r, 2])) #trapezoid rule; this calculates area of each chunk
            CURVE <- c(CURVE, curve)
          } #write down each chunk
          
          area <- sum(CURVE) #take the sum
          output <- c(i, j, p, m, n, area)
          OUT <- rbind(OUT, output)
          
        }
      }
    }
  }
}




colnames(OUT) <- c("Year", "Lake", "Parasite.Species", "Focal.Host", "Other.Host", "AUC.prev")
# write.table(OUT, "auc_14to16_prev.txt", col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE) 

# tidy up the data frame
auc_shared <- as.data.frame(OUT)
row.names(auc_shared) <- c()
auc_shared$AUC.prev <- as.numeric(as.character(auc_shared$AUC.prev))

lakeyear_auc_shared_prev <- auc_shared

write.csv(lakeyear_auc_shared_prev, "lakeyear_auc_shared_prev_geo.csv")

#===================================================================================
# SHARED AUC PREVALENCE*DENSITY
#===================================================================================
# general code ---------------------------------------------------------------------
clean_auc <- read.csv("Clean-Data-2014-2016.csv", header = TRUE)

# rename columns; easier to "gather" this way
colnames(clean_auc)[which(names(clean_auc) == "pasteuria.prev")] <- "pasteuria"
colnames(clean_auc)[which(names(clean_auc) == "metsch.prev")] <- "metsch"
colnames(clean_auc)[which(names(clean_auc) == "scarlet.prev")] <- "scarlet"
colnames(clean_auc)[which(names(clean_auc) == "brood.prev")] <- "brood"
colnames(clean_auc)[which(names(clean_auc) == "larssonia.prev")] <- "larssonia"
colnames(clean_auc)[which(names(clean_auc) == "spider.prev")] <- "spider"
colnames(clean_auc)[which(names(clean_auc) == "gurleya.prev")] <- "gurleya"

# split data depending on total number of hosts counted in parasite sample
total_less <- clean_auc %>%
  filter(Total < 20)

# if less than 20 total, set prevalence to 0
total_less$pasteuria <- 0
total_less$metsch <- 0
total_less$scarlet <- 0
total_less$brood <- 0
total_less$larssonia <- 0
total_less$spider <- 0
total_less$gurleya <- 0

# don't do anything if total is greater/equal to 20
total_greater <- clean_auc %>%
  filter(Total >= 20)

# join all data back together
loop_data <- full_join(total_less, total_greater)

# gather into tall format; select appropriate columns
loop_data <- loop_data %>%
  gather(Parasite.Species, Prevalence, pasteuria:gurleya) %>%
  select(Host.Species, Year, Lake, Parasite.Species, Julian.Day, Prevalence, Density) %>%
  mutate(Prev.Density = Prevalence*Density)


loop_data <- loop_data %>%
  select(Host.Species, Year, Lake, Parasite.Species, Julian.Day, Prev.Density)

# set year and julian day as correct class
loop_data$Year <- as.factor(loop_data$Year)
loop_data$Julian.Day <- as.integer(loop_data$Julian.Day)

# order everything by julian day
loop_data <- loop_data[order(loop_data$Julian.Day), ] #make sure julian days are in order

# make sure it's in data frame format
loop_data <- as.data.frame(loop_data)


# Shared AUC loop -------------------------------------------------------------------------
# calculates the overlapping prevalence area (same parasite, different hosts)

YEAR<-unique(loop_data$Year)

OUT <- NULL

for (i in YEAR) {
  thisyear <- loop_data[loop_data$Year == i, ] # look at 1 year at a time
  LAKE <- unique(thisyear$Lake)
  
  for (j in LAKE) {
    thislake <- thisyear[thisyear$Lake == j, ] # within a year, look at 1 lake at a time
    PARA <- unique(thislake$Parasite.Species)
    
    for (p in PARA) {
      thisparasite <- thislake[thislake$Parasite.Species == p, ] # within a lake, look at 1 parasite at a time
      SPECIES <- unique(thisparasite$Host.Species)
      
      for (m in SPECIES) { #where m = dent, retro, cerio, etc
        thisspecies <- thisparasite[thisparasite$Host.Species == m, ]
        stuffone <- select(thisspecies, Julian.Day, Prev.Density) # focal host data
        
        for (n in SPECIES) {
          otherspecies <- thisparasite[thisparasite$Host.Species == n, ]
          stufftwo <- select(otherspecies, Julian.Day, Prev.Density) # other host data
          both <- full_join(stuffone, stufftwo, by = "Julian.Day") # use whatever merging function makes sense
          
          
          both$Prev.Density.x[both$Prev.Density.x == 0] <- NA  
          both$Prev.Density.y[both$Prev.Density.y == 0] <- NA  
          avg_epi <- rowMeans(both[, c("Prev.Density.x", "Prev.Density.y")], na.rm = FALSE)
          both$avg_epi <- avg_epi   
          both[is.na(both)] <- 0
          avgboth <- both
          # avgboth will have columns: Julian.Day, Prev.Density.x, Prev.Density.y, avg_epi
          TOTS <- NULL #make an empty matrix to include sequential Prev.Densitys. This will be used to calculate integrated areas.
          
          # now calculate the auc using the avgboth data
          for (q in 1:length(avgboth$Prev.Density.x)) { #CHECK THAT THIS PART MAKES SENSE
            total <- avgboth[q, 4] # average Prev.Density Prev.Density
            jul <- avgboth[q, 1] # julian day
            tots <- c(q, jul, total) # new row in dataframe with q, julian day, and total
            TOTS <- rbind(TOTS, tots) # bind all rows together (different one for each lake/year/parasite)
            colnames(TOTS) <- c("q", "date", "total")
          }
          
          rownames(TOTS) <- NULL
          TOTS <- as.data.frame(TOTS)
          CURVE <- NULL #make an empty matrix for calulating each chunk of the integrated area by trapezoid rule.
          
          for (r in 1:(length(TOTS$total) - 1)) {
            curve <- 0.5*((TOTS[r, 3] + TOTS[r+1, 3])*(TOTS[r+1, 2] - TOTS[r, 2])) #trapezoid rule; this calculates area of each chunk
            CURVE <- c(CURVE, curve)
          } #write down each chunk
          
          area <- sum(CURVE) #take the sum
          output <- c(i, j, p, m, n, area)
          OUT <- rbind(OUT, output)
          
        }
      }
    }
  }
}




colnames(OUT) <- c("Year", "Lake", "Parasite.Species", "Focal.Host", "Other.Host", "AUC.prev")
# write.table(OUT, "auc_14to16_prev.txt", col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE) 

# tidy up the data frame
auc_shared <- as.data.frame(OUT)
row.names(auc_shared) <- c()
auc_shared$AUC.prev <- as.numeric(as.character(auc_shared$AUC.prev))

lakeyear_auc_shared_prevdens <- auc_shared

# Shared AUC loop (geometric mean) -------------------------------------------------------------------------
# calculates the overlapping prevalence area (same parasite, different hosts)

YEAR<-unique(loop_data$Year)

OUT <- NULL

for (i in YEAR) {
  thisyear <- loop_data[loop_data$Year == i, ] # look at 1 year at a time
  LAKE <- unique(thisyear$Lake)
  
  for (j in LAKE) {
    thislake <- thisyear[thisyear$Lake == j, ] # within a year, look at 1 lake at a time
    PARA <- unique(thislake$Parasite.Species)
    
    for (p in PARA) {
      thisparasite <- thislake[thislake$Parasite.Species == p, ] # within a lake, look at 1 parasite at a time
      SPECIES <- unique(thisparasite$Host.Species)
      
      for (m in SPECIES) { #where m = dent, retro, cerio, etc
        thisspecies <- thisparasite[thisparasite$Host.Species == m, ]
        stuffone <- select(thisspecies, Julian.Day, Prev.Density) # focal host data
        
        for (n in SPECIES) {
          otherspecies <- thisparasite[thisparasite$Host.Species == n, ]
          stufftwo <- select(otherspecies, Julian.Day, Prev.Density) # other host data
          both <- full_join(stuffone, stufftwo, by = "Julian.Day") # use whatever merging function makes sense
          
          
          both$Prev.Density.x[both$Prev.Density.x == 0] <- NA  
          both$Prev.Density.y[both$Prev.Density.y == 0] <- NA  
          avg_epi <- apply(both[, c("Prev.Density.x", "Prev.Density.y")], 1, function(x) (prod(x[x!=0]))^(1/sum(x!=0)))
          both$avg_epi <- avg_epi   
          both[is.na(both)] <- 0
          avgboth <- both
          # avgboth will have columns: Julian.Day, Prev.Density.x, Prev.Density.y, avg_epi
          TOTS <- NULL #make an empty matrix to include sequential Prev.Densitys. This will be used to calculate integrated areas.
          
          # now calculate the auc using the avgboth data
          for (q in 1:length(avgboth$Prev.Density.x)) { #CHECK THAT THIS PART MAKES SENSE
            total <- avgboth[q, 4] # average Prev.Density Prev.Density
            jul <- avgboth[q, 1] # julian day
            tots <- c(q, jul, total) # new row in dataframe with q, julian day, and total
            TOTS <- rbind(TOTS, tots) # bind all rows together (different one for each lake/year/parasite)
            colnames(TOTS) <- c("q", "date", "total")
          }
          
          rownames(TOTS) <- NULL
          TOTS <- as.data.frame(TOTS)
          CURVE <- NULL #make an empty matrix for calulating each chunk of the integrated area by trapezoid rule.
          
          for (r in 1:(length(TOTS$total) - 1)) {
            curve <- 0.5*((TOTS[r, 3] + TOTS[r+1, 3])*(TOTS[r+1, 2] - TOTS[r, 2])) #trapezoid rule; this calculates area of each chunk
            CURVE <- c(CURVE, curve)
          } #write down each chunk
          
          area <- sum(CURVE) #take the sum
          output <- c(i, j, p, m, n, area)
          OUT <- rbind(OUT, output)
          
        }
      }
    }
  }
}




colnames(OUT) <- c("Year", "Lake", "Parasite.Species", "Focal.Host", "Other.Host", "AUC.prev")
# write.table(OUT, "auc_14to16_prev.txt", col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE) 

# tidy up the data frame
auc_shared <- as.data.frame(OUT)
row.names(auc_shared) <- c()
auc_shared$AUC.prev <- as.numeric(as.character(auc_shared$AUC.prev))

lakeyear_auc_shared_prevdensgeo <- auc_shared

write.csv(lakeyear_auc_shared_prevdensgeo, "lakeyear_auc_shared_prevdens_geo.csv")


# centrality-lake-figs

lake_auc_density <- lakeyear_auc_density %>%
  group_by(Lake, Host.Species) %>%
  summarise(AUC.density = sum(AUC.density)) %>%
  mutate(Log.AUC.density = log(AUC.density + 1))


# EPIDEMIC OVERLAP--PREVALENCE ===========================================================

# arrange data for make_mat_list
auc_shared <- lakeyear_auc_shared_prev %>%
  group_by(Lake, Parasite.Species, Focal.Host, Other.Host) %>%
  summarise(Edge.Weight = sum(AUC.prev))
auc_shared <- as.data.frame(auc_shared)

# MAD is adding these next bits to convert from character to factor because it's not working with them as characters
auc_shared$Parasite.Species <- factor(auc_shared$Parasite.Species)
auc_shared$Lake <- factor(auc_shared$Lake)
auc_shared$Focal.Host <- factor(auc_shared$Focal.Host)

parasite_vec <- c(levels(auc_shared$Parasite.Species)) # set up vectors of names for loop
host_vec <- c(levels(auc_shared$Focal.Host))
lake_vec <- c(levels(auc_shared$Lake))
lakes <- length(lake_vec)
paras <- length(parasite_vec)
paralake <- lakes*paras

#names(auc_shared)[names(auc_shared) == 'AUC.prev'] <- 'Edge.Weight'

# make list of host adjacency matrics
auc_shared_mat <- make_mat_list(auc_shared, parasite_vec, lake_vec) 
names(auc_shared_mat) <- parasite_vec

# calculate centrality
raw_centrality_data <- data.frame(matrix(rep(NA, 1*13), nrow = 1)) #outside of loop, set up empty matrix for data 
names(raw_centrality_data) <- c(host_vec, "Lake", "Parasite.Species",  "Metric", "Vertices", "Edges")

centrality_bylake <- calc_centrality_bylake(auc_shared_mat, raw_centrality_data, parasite_vec, lake_vec)

# clean/organize centrality data
rank_lake_central <- organize_centrality(centrality_bylake)

# calculate and plot graph density
raw_graphdensity_data <- data.frame(matrix(rep(NA, paralake*5), nrow = paralake))
names(raw_graphdensity_data) <- c("Parasite.Species", "Lake","Vertices", "Edges", "Graph.Density")

graph_density_data <- calc_graph_density(auc_shared_mat, raw_graphdensity_data, parasite_vec, lake_vec, lake_auc_density)
graph_density_plots <- plot_graph_density(graph_density_data)
names(graph_density_plots) <- c("all_data", "nonzero_data")

# clone a separate df so we can compare across edge weight types
prevalence_rank <- rank_lake_central

# EPIDEMIC OVERLAP--INFECTED HOST DENSITY ===========================================================

#arrange data for make_mat_list (don't log-transform version)
auc_shared <- lakeyear_auc_shared_prevdens %>%
  group_by(Lake, Parasite.Species, Focal.Host, Other.Host) %>%
  summarise(Edge.Weight = sum(AUC.prev)) 
auc_shared$Lake <- as.factor(auc_shared$Lake)
auc_shared <- as.data.frame(auc_shared)

# MAD is adding these next bits to convert from character to factor because it's not working with them as characters
auc_shared$Parasite.Species <- factor(auc_shared$Parasite.Species)
auc_shared$Lake <- factor(auc_shared$Lake)
auc_shared$Focal.Host <- factor(auc_shared$Focal.Host)

parasite_vec <- c(levels(auc_shared$Parasite.Species)) # set up vectors of names for loop
host_vec <- c(levels(auc_shared$Focal.Host))
lake_vec <- c(levels(auc_shared$Lake))

# make list of host adjacency matrics
auc_shared_mat <- make_mat_list(auc_shared, parasite_vec, lake_vec) 
names(auc_shared_mat) <- parasite_vec

# calculate centrality
raw_centrality_data <- data.frame(matrix(rep(NA, 1*13), nrow = 1)) #outside of loop, set up empty matrix for data 
names(raw_centrality_data) <- c(host_vec, "Lake", "Parasite.Species",  "Metric", "Vertices", "Edges")

centrality_bylake <- calc_centrality_bylake(auc_shared_mat, raw_centrality_data, parasite_vec, lake_vec)

# clean/organize centrality data
rank_lake_central <- organize_centrality(centrality_bylake)

# clone a separate df so we can compare across edge weight types
infecteds_rank <- rank_lake_central

# COMPARING ACROSS METRICS ===========================================================

# dfs with centrality scores for each edge weight method:
# prevalence_rank 
# infecteds_rank

# make a new factor column for distinguishing data used to create edge weights
prevalence_data <- prevalence_rank %>% 
  mutate(Weight.Method = factor("prevalence",
                                levels = c("prevalence", "infecteds")))
infecteds_data <- infecteds_rank %>% 
  mutate(Weight.Method = factor("infecteds",
                                levels = c("prevalence", "infecteds")))

lake_diff_edges <- bind_rows(prevalence_data, infecteds_data)
common_parasite_vec <- c("brood", "scarlet", "pasteuria")

# GRAPH DENSITY ======================================

# calculate graph density stats
nonzero_graph_density <- graph_density_data %>% filter(Graph.Density > 0) # filter data
nonzero_graph_density$Parasite.Species <- as.factor(nonzero_graph_density$Parasite.Species)
nonzero_graph_density$Lake <- as.factor(nonzero_graph_density$Lake)
nonzero_graph_density$Vertices <- as.integer(nonzero_graph_density$Vertices)
nonzero_graph_density$Edges <- as.integer(nonzero_graph_density$Edges)
nonzero_graph_density$Graph.Density <- as.numeric(nonzero_graph_density$Graph.Density)

raw_analyses_parasites <- data.frame(matrix(rep(NA, 21*4), nrow = 21)) 
names(raw_analyses_parasites) <- c( "KW.p", "KW.chisq", "Comparison", "Dunn.p")

graph_density_stats <- compare_parasites(nonzero_graph_density, raw_analyses_parasites)

graph_density_sig <- graph_density_stats %>%
  filter(Dunn.p < 0.05)

gd_kwtest <- min(graph_density_stats$KW.p)

kw_result <- kruskal.test(Graph.Density ~ Parasite.Species, data = nonzero_graph_density)
dunn_result <- dunn.test(nonzero_graph_density$Graph.Density, nonzero_graph_density$Parasite.Species,
                         kw = TRUE, method = "bonferroni", altp = TRUE)
coniman_result <- conover.test(nonzero_graph_density$Graph.Density, nonzero_graph_density$Parasite.Species, kw = TRUE, 
                               method = "bonferroni", altp = TRUE)

# PLOT 1: graph density ================================================================

# Graph density with conover-iman post-hoc
max_gdens <- nonzero_graph_density %>%
  group_by(Parasite.Species) %>%
  summarise(max = max(Graph.Density)) %>%
  arrange(max)
ytext_vec <- c(max_gdens$max + 0.125)

xtext_vec <- seq(1, 7, by = 1)
mtext_vec <- c("bc", "c", "bc", "bc", "ab", "ab", "a")
ytext_vec <- c(max_gdens$max + 0.125)


nonzero_graph_density$Parasite.Species <- str_replace(nonzero_graph_density$Parasite.Species, "scarlet", "spiro")

gdensity_plot <- ggplot(data = nonzero_graph_density, 
                        aes(x = reorder(Parasite.Species, Graph.Density, FUN = median), 
                            y = Graph.Density)) +
  geom_boxplot(fill = "#f0f0f0", outlier.shape = NA) +
  geom_jitter(data = nonzero_graph_density, aes(x = reorder(Parasite.Species, Graph.Density, FUN = median),
                                                y = Graph.Density),
              width = 0.1, height = 0, size = 2.5, alpha = 0.4) +
  xlab("") +
  scale_x_discrete(labels=c("G.vavrai", "L.obtusa", "Spider","M.bicuspidata","P.ramosa","S.cienkowskii","B.paedophthorum"))+
  ylab("Graph density \n (total edges / total possible edges)") +
  ylim(0, 1) +
  ggtitle("") +
  #annotate("text", x = 0.75, y = 0.75, label = paste("Kruskal-Wallis, p < 0.0001")) +
  #annotate("text", x = 1, y = 0.75, label = "Pairwise test: Conover-Iman") +
  annotate("text", x = xtext_vec, y = ytext_vec, label = c(mtext_vec)) +
  coord_flip()
gdensity_plot

pdf(file = paste("fig_graphdensity.pdf"), width = 4.25, height = 4, useDingbats=F)
gdensity_plot
dev.off()

# get means (for references in paper)
mean_gdensity <- nonzero_graph_density %>% 
  group_by(Parasite.Species) %>%
  summarise(mn = mean(Graph.Density))



# PLOT 2: host centrality ================================================================


# (B.1) Eigenvector centrality--Pasteuria prevalence
xtext_vec <- c(seq(1, 5, by = 1))
mtext_vec <- c("b", "b", "ab", "a", "a")
ytext_vec <- c(rep(1.12, 5))

length(ytext_vec)

subset_eigen <- lake_diff_edges %>%
  filter(Metric == "eigenvector") %>%
  filter(Parasite.Species == "pasteuria") %>%
  filter(Weight.Method == "prevalence")

eigen_past_prev_plot <- ggplot() +
  geom_boxplot(data = subset_eigen, aes(x = reorder(Host.Species, Centrality, FUN = median), y = Centrality),
               fill = "#377eb8", outlier.shape = NA, alpha = 0.4) +
  geom_jitter(data = subset_eigen, aes(x = reorder(Host.Species, Centrality, FUN = median), 
                                       y = Centrality), color = "black",
              width = .1, height = 0, size = 2.5, alpha = 0.4) +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.0),
                     labels = c("0", "0.25", "0.5", "0.75", "1"),
                     limits = c(0, 1.125)) +
  xlab("") +
  ylab("") +
  annotate("text", x = xtext_vec, y = ytext_vec, label = c(mtext_vec)) +
  coord_flip() +
  theme_bw()
eigen_past_prev_plot

# (B.2) Eigenvector centrality--Pasteuria infecteds
xtext_vec <- c(seq(1, 5, by = 1))
mtext_vec <- c("c", "bc", "bc", "ab", "ab")
ytext_vec <- c(rep(1.12, 5))

subset_eigen <- lake_diff_edges %>%
  filter(Metric == "eigenvector") %>%
  filter(Parasite.Species == "pasteuria") %>%
  filter(Weight.Method == "infecteds")

eigen_past_inf_plot <- ggplot() +
  geom_boxplot(data = subset_eigen, aes(x = reorder(Host.Species, Centrality, FUN = median), y = Centrality),
               fill = "#377eb8", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(data = subset_eigen, aes(x = reorder(Host.Species, Centrality, FUN = median), 
                                       y = Centrality), color = "black",
              width = .1, height = 0, size = 2.5, alpha = 0.4) +
  ylab("") +
  xlab("") +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.0),
                     labels = c("0", "0.25", "0.5", "0.75", "1"),
                     limits = c(0, 1.125)) +
  annotate("text", x = xtext_vec, y = ytext_vec, label = c(mtext_vec)) +
  coord_flip() +
  theme_bw()
eigen_past_inf_plot


# (C.1) Eigenvector centrality--Brood prevalence
subset_eigen <- lake_diff_edges %>%
  filter(Metric == "eigenvector") %>%
  filter(Parasite.Species == "brood") %>%
  filter(Weight.Method == "prevalence")

eigen_brood_prev_plot <- ggplot() +
  geom_boxplot(data = subset_eigen, aes(x = reorder(Host.Species, Centrality, FUN = median), y = Centrality),
               fill = "#4daf4a", outlier.shape = NA, alpha = 0.4) +
  geom_jitter(data = subset_eigen, aes(x = reorder(Host.Species, Centrality, FUN = median), 
                                       y = Centrality), color = "black",
              width = .1, height = 0, size = 2.5, alpha = 0.4) +
  ylab("") +
  xlab("") +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.0),
                     labels = c("0", "0.25", "0.5", "0.75", "1"),
                     limits = c(0, 1.125)) +
  coord_flip() +
  theme_bw()
eigen_brood_prev_plot

# (C.2) Eigenvector centrality--Brood infecteds
xtext_vec <- c(0.75, 1.75, 2.75, 3.75, 5, 6)
mtext_vec <- c("b", "b", "ab", "ab", "ab", "a")
ytext_vec <- c(rep(1.12, 6))

subset_eigen <- lake_diff_edges %>%
  filter(Metric == "eigenvector") %>%
  filter(Parasite.Species == "brood") %>%
  filter(Weight.Method == "infecteds")

eigen_brood_inf_plot <- ggplot() +
  geom_boxplot(data = subset_eigen, aes(x = reorder(Host.Species, Centrality, FUN = median), y = Centrality),
               fill = "#4daf4a", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(data = subset_eigen, aes(x = reorder(Host.Species, Centrality, FUN = median), 
                                       y = Centrality), color = "black",
              width = .1, height = 0, size = 2.5, alpha = 0.4) +
  ylab("") +
  xlab("") +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.0),
                     labels = c("0", "0.25", "0.5", "0.75", "1"),
                     limits = c(0, 1.125)) +
  annotate("text", x = xtext_vec, y = ytext_vec, label = c(mtext_vec)) +
  coord_flip() +
  theme_bw()
eigen_brood_inf_plot

# (D.1) Eigenvector centrality--Scarlet prevalence
subset_eigen <- lake_diff_edges %>%
  filter(Metric == "eigenvector") %>%
  filter(Parasite.Species == "scarlet") %>%
  filter(Weight.Method == "prevalence")

eigen_scar_prev_plot <- ggplot() +
  geom_boxplot(data = subset_eigen, aes(x = reorder(Host.Species, Centrality, FUN = median), y = Centrality),
               fill = "#e41a1c", outlier.shape = NA, alpha = 0.4) +
  geom_jitter(data = subset_eigen, aes(x = reorder(Host.Species, Centrality, FUN = median), 
                                       y = Centrality), color = "black",
              width = .1, height = 0, size = 2.5, alpha = 0.4) +
  ylab("Eigenvector centrality") +
  xlab("") +
  #ggtitle("Spiro \noverlapping prevalence") +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.0),
                     labels = c("0", "0.25", "0.5", "0.75", "1"),
                     limits = c(0, 1.125)) +
  coord_flip() +
  theme_bw()
eigen_scar_prev_plot

# (D.2) Eigenvector centrality--Scarlet infecteds
subset_eigen <- lake_diff_edges %>%
  filter(Metric == "eigenvector") %>%
  filter(Parasite.Species == "scarlet") %>%
  filter(Weight.Method == "infecteds")

eigen_scar_inf_plot <- ggplot() +
  geom_boxplot(data = subset_eigen, aes(x = reorder(Host.Species, Centrality, FUN = median), y = Centrality),
               fill = "#e41a1c", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(data = subset_eigen, aes(x = reorder(Host.Species, Centrality, FUN = median), 
                                       y = Centrality), color = "black",
              width = .1, height = 0, size = 2.5, alpha = 0.4) +
  ylab("Eigenvector centrality") +
  xlab("") +
  #ggtitle("Spiro \noverlapping infected host density") +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.0),
                     labels = c("0", "0.25", "0.5", "0.75", "1"),
                     limits = c(0, 1.125)) +
  # ggtitle("Scarlet Infecteds") +
  # annotate("text", x = 6.5, y = 0.25, label = paste("Kruskal-Wallis, p = 0.5")) +
  coord_flip() + 
  theme_bw()
eigen_scar_inf_plot

# plot all plots together
centr_plot <- plot_grid(eigen_past_prev_plot, eigen_past_inf_plot,
                        eigen_brood_prev_plot, eigen_brood_inf_plot,
                        eigen_scar_prev_plot, eigen_scar_inf_plot,
                        ncol = 2, labels = "AUTO")

# add labels for parasites to rows
centr_labs_right <- ggdraw() +
  draw_label("P.ramosa", x = 0.2, y = 0.88, fontface = "italic", angle = 270) +
  draw_label("B.paedophthorum", x = 0.2, y = 0.55, fontface = "italic", angle = 270) +
  draw_label("S.cienkowskii", x = 0.2, y = 0.22, fontface = "italic", angle = 270)

# add labels for metric to columns
centr_labs_top <- ggdraw() +
  draw_label("Prevalence", x = 0.25, y = 0.2, fontface = "italic", angle = 0) +
  draw_label("Infected host density", x = 0.72, y = 0.2, fontface = "italic", angle = 0) 

# combine labels independently
centr_right <- plot_grid(centr_plot, centr_labs_right, rel_widths = c(1, 0.1))  
centr_top <- plot_grid(centr_labs_top, centr_right, ncol = 1, rel_heights  = c(0.1, 1))  

pdf(file = "fig_centrality_avg_epi_overlap.pdf", width = 7, height = 8, useDingbats=F)
centr_top # save
dev.off()


# STATS: centrality ================================================================


# Kruskal-Wallis test & Conover-Iman post hoc
edge_method_vec <- c("prevalence", "infecteds")

raw_analyses_hosts <- data.frame(matrix(rep(NA, 120*6), nrow = 120)) 
names(raw_analyses_hosts) <- c("Parasite.Species", "Weight.Method", "KW.p", "KW.chisq", 
                               "Comparison", "Conover.p")
centrality_stats_conover <- compare_hosts_con(lake_diff_edges, common_parasite_vec, edge_method_vec, raw_analyses_hosts)

centrality_kw_sig <- centrality_stats_conover %>%
  filter(KW.p < 0.05)

centrality_conover_sig <- centrality_stats_conover %>%
  filter(Conover.p < 0.05)

write.csv(centrality_conover_sig, "centrality_conover_sig_avg_epi.csv")


# PLOTS: AUC single epidemic size ===========================================================

epidemic_size <- lakeyear_auc_shared_prev %>%
  filter(Focal.Host == Other.Host) %>%
  filter(Focal.Host != "mendotae") %>%
  filter(AUC.prev > 0) %>%
  group_by(Parasite.Species, Focal.Host) %>% 
  mutate(AUC.log = log(AUC.prev + 1)) %>%
  select(Parasite.Species, Focal.Host, AUC.prev, AUC.log)
epidemic_size <- as.data.frame(epidemic_size)

brood_prev_epi <- epidemic_size %>% filter(Parasite.Species == "brood")
scarlet_prev_epi <- epidemic_size %>% filter(Parasite.Species == "scarlet")
past_prev_epi <- epidemic_size %>% filter(Parasite.Species == "pasteuria")

epi_past_prev_plot <- ggplot() +
  geom_boxplot(data = past_prev_epi, aes(x = reorder(Focal.Host, AUC.prev, FUN = median), y = AUC.prev),
               fill = "#377eb8", outlier.shape = NA, alpha = 0.4) +
  geom_jitter(data = past_prev_epi, aes(x = reorder(Focal.Host, AUC.prev, FUN = median), 
                                        y = AUC.prev),
              width = .1, height = 0, size = 2, alpha = 0.4) +
  ylab("") +
  xlab("") +
  # ylim(0, 18) +
  coord_flip() +
  theme_bw()
epi_past_prev_plot

epi_brood_prev_plot <- ggplot() +
  geom_boxplot(data = brood_prev_epi, aes(x = reorder(Focal.Host, AUC.prev, FUN = median), y = AUC.prev),
               fill = "#4daf4a", outlier.shape = NA, alpha = 0.4) +
  geom_jitter(data = brood_prev_epi, aes(x = reorder(Focal.Host, AUC.prev, FUN = median), 
                                         y = AUC.prev),
              width = .1, height = 0, size = 2, alpha = 0.4) +
  ylab("") +
  xlab("") +
  # ylim(0, 18) +
  coord_flip() +
  theme_bw()
epi_brood_prev_plot

epi_scar_prev_plot <- ggplot() +
  geom_boxplot(data = scarlet_prev_epi, aes(x = reorder(Focal.Host, AUC.prev, FUN = median), y = AUC.prev),
               fill = "#e41a1c", outlier.shape = NA, alpha = 0.4) +
  geom_jitter(data = scarlet_prev_epi, aes(x = reorder(Focal.Host, AUC.prev, FUN = median), 
                                           y = AUC.prev),
              width = .1, height = 0, size = 2, alpha = 0.4) +
  ylab("Epidemic size (AUC)") +
  xlab("") +
  # ylim(0, 18) +
  coord_flip() +
  theme_bw()
epi_scar_prev_plot

# stats for single epidemic AUC of prevalence
scar_epid_kw <- kruskal.test(AUC.prev ~ Focal.Host, data = scarlet_prev_epi)
brood_epid_kw <- kruskal.test(AUC.prev ~ Focal.Host, data = brood_prev_epi)
past_epid_kw <- kruskal.test(AUC.prev ~ Focal.Host, data = past_prev_epi)

# PLOTS: AUC single epidemic infected host density ============================================

epidemic_size_inf <- lakeyear_auc_shared_prevdens %>%
  filter(Focal.Host == Other.Host) %>%
  filter(Focal.Host != "mendotae") %>%
  filter(AUC.prev > 0) %>%
  group_by(Parasite.Species, Focal.Host) %>% 
  mutate(AUC.log = log(AUC.prev + 1)) %>%
  select(Parasite.Species, Focal.Host, AUC.prev, AUC.log)
epidemic_size_inf <- as.data.frame(epidemic_size_inf)

brood_prev_epi_inf <- epidemic_size_inf %>% filter(Parasite.Species == "brood")
scarlet_prev_epi_inf <- epidemic_size_inf %>% filter(Parasite.Species == "scarlet")
past_prev_epi_inf <- epidemic_size_inf %>% filter(Parasite.Species == "pasteuria")

common_epis <-epidemic_size_inf %>% filter(Parasite.Species == "brood"|
                                             Parasite.Species == "scarlet"|
                                             Parasite.Species == "pasteuria")
min(common_epis$AUC.log)
max(common_epis$AUC.log)


xtext_vec <- c(0.75, 1.75, 2.75, 3.75, 5, 6)
mtext_vec <- c("b", "ab", "ab", "ab", "a", "a")
ytext_vec <- c(rep(17, 6))

epi_past_inf_plot <- ggplot() +
  geom_boxplot(data = past_prev_epi_inf, aes(x = reorder(Focal.Host, AUC.log, FUN = median), y = AUC.log),
               fill = "#377eb8", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(data = past_prev_epi_inf, aes(x = reorder(Focal.Host, AUC.log, FUN = median), 
                                            y = AUC.log),
              width = .1, height = 0, size = 2, alpha = 0.4) +
  ylab("") +
  xlab("") +
  ylim(6, 17.5) +
  annotate("text", x = xtext_vec, y = ytext_vec, label = c(mtext_vec)) +
  coord_flip() +
  theme_bw()
epi_past_inf_plot


xtext_vec <- c(0.75, 1.75, 2.75, 3.75, 5, 6)
mtext_vec <- c("c", "bc", "abc", "ab", "a", "ab")
ytext_vec <- c(rep(17, 6))

epi_brood_inf_plot <- ggplot() +
  geom_boxplot(data = brood_prev_epi_inf, aes(x = reorder(Focal.Host, AUC.log, FUN = median), y = AUC.log),
               fill = "#4daf4a", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(data = brood_prev_epi_inf, aes(x = reorder(Focal.Host, AUC.log, FUN = median), 
                                             y = AUC.log),
              width = .1, height = 0, size = 2, alpha = 0.4) +
  ylab("") +
  xlab("") +
  ylim(6, 17.5) +
  annotate("text", x = xtext_vec, y = ytext_vec, label = c(mtext_vec)) +
  coord_flip() +
  theme_bw()
epi_brood_inf_plot

epi_scar_inf_plot <- ggplot() +
  geom_boxplot(data = scarlet_prev_epi_inf, aes(x = reorder(Focal.Host, AUC.log, FUN = median), y = AUC.log),
               fill = "#e41a1c", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(data = scarlet_prev_epi_inf, aes(x = reorder(Focal.Host, AUC.log, FUN = median), 
                                               y = AUC.log),
              width = .1, height = 0, size = 2, alpha = 0.4) +
  ylab("Epidemic size (log AUC)") +
  xlab("") +
  ylim(6, 17.5) +
  coord_flip() +
  theme_bw()
epi_scar_inf_plot


auc_plot <- plot_grid(epi_past_prev_plot, epi_past_inf_plot,
                      epi_brood_prev_plot, epi_brood_inf_plot,
                      epi_scar_prev_plot, epi_scar_inf_plot,
                      ncol = 2, labels = "AUTO")

# add labels for parasites to rows
auc_labs_right <- ggdraw() +
  draw_label("P.ramosa", x = 0.2, y = 0.88, fontface = "italic", angle = 270) +
  draw_label("B.paedophthorum", x = 0.2, y = 0.55, fontface = "italic", angle = 270) +
  draw_label("S.cienkowskii", x = 0.2, y = 0.22, fontface = "italic", angle = 270)

# add labels for metric to columns
auc_labs_top <- ggdraw() +
  draw_label("Prevalence", x = 0.25, y = 0.2, fontface = "italic", angle = 0) +
  draw_label("Infected host density", x = 0.72, y = 0.2, fontface = "italic", angle = 0) 

# combine labels independently
auc_right <- plot_grid(auc_plot, auc_labs_right, rel_widths = c(1, 0.1))  
auc_top <- plot_grid(auc_labs_top, auc_right, ncol = 1, rel_heights  = c(0.1, 1)) 

pdf(file = "fig_AUC.pdf", width = 6.75, height = 8, useDingbats=F)
auc_top
dev.off()



# stats for single epidemic AUC of prevalence
scar_epi_inf_kw <- kruskal.test(AUC.log ~ Focal.Host, data = scarlet_prev_epi_inf)
brood_epi_inf_kw <- kruskal.test(AUC.log ~ Focal.Host, data = brood_prev_epi_inf)
past_epi_inf_kw <- kruskal.test(AUC.log ~ Focal.Host, data = past_prev_epi_inf)

brood_con_i <- conover.test(brood_prev_epi_inf$AUC.log, brood_prev_epi_inf$Focal.Host, kw = TRUE, 
                            method = "bonferroni", altp = TRUE)
past_con_i <- conover.test(past_prev_epi_inf$AUC.log, past_prev_epi_inf$Focal.Host, kw = TRUE, 
                           method = "bonferroni", altp = TRUE)



# host total density ==============================================
str(lake_auc_density)
lake_auc_density <- as.data.frame(lake_auc_density)

# only include lakes where host species was present
auc_nonzero_density <- lake_auc_density %>% 
  filter(AUC.density > 0) %>%
  filter(Host.Species != "mendotae")

total_host_auc_plot <- ggplot() +
  geom_boxplot(data = auc_nonzero_density, aes(x = reorder(Host.Species, AUC.density, FUN = median), y = AUC.density/100000000),
               fill = "black", outlier.shape = NA, alpha = 0.5) +
  geom_jitter(data = auc_nonzero_density, aes(x = reorder(Host.Species, AUC.density, FUN = median), 
                                               y = AUC.density/100000000),
              width = .1, height = 0, size = 2, alpha = 0.4) +
  ylab("Host density (AUC * 10^8)") +
  xlab("") +
  coord_flip() 
total_host_auc_plot

logtotal_host_auc_plot <- ggplot() +
  geom_boxplot(data = auc_nonzero_density, aes(x = reorder(Host.Species, Log.AUC.density, FUN = median), y = Log.AUC.density),
               fill = "black", outlier.shape = NA, alpha = 0.5) +
  geom_jitter(data = auc_nonzero_density, aes(x = reorder(Host.Species, Log.AUC.density, FUN = median), 
                                           y = Log.AUC.density),
              width = .1, height = 0, size = 2, alpha = 0.4) +
  ylab("Host density (log AUC)") +
  xlab("") +
  coord_flip() 
logtotal_host_auc_plot


host_auc_tog <- plot_grid(total_host_auc_plot, logtotal_host_auc_plot, ncol = 2) 

pdf(file = "fig_host_density_AUC.pdf", width = 8, height = 4, useDingbats=F)
host_auc_tog
dev.off()


# KW test
host_auc_kw <- kruskal.test(Log.AUC.density ~ Host.Species, data = auc_nonzero_density)
# Con-Iman post hoc
host_auc_con_i <- conover.test(auc_nonzero_density$Log.AUC.density, auc_nonzero_density$Host.Species, kw = TRUE, 
                            method = "bonferroni", altp = TRUE)


# messing around space ==============================================

str(prevalence_rank)
str(lake_auc_density)

eigen_prev <- prevalence_rank %>% filter(Metric == "eigenvector")

cent_v_density <- full_join(eigen_prev, lake_auc_density)

cent_v_density <- cent_v_density %>% na.omit()
cent_v_density$Host.Species <- as.factor(cent_v_density$Host.Species)

str(cent_v_density)

p1 <- ggplot(data = cent_v_density %>% filter(Parasite.Species == "brood"|
                                                Parasite.Species == "scarlet"|
                                                Parasite.Species == "pasteuria"),
             aes(x = AUC.density, y = Centrality, color = Host.Species)) +
  geom_point() +
  facet_grid(Host.Species~Parasite.Species, scales = "free")
p1
















# data with main parasite species and AUCs in different hosts
host_auc_prev <- lakeyear_auc_shared_prev %>%
  filter(Focal.Host == Other.Host) %>%
  group_by(Lake, Focal.Host, Parasite.Species) %>%
  summarise(AUC.single = sum(AUC.prev)) %>%
  filter(Parasite.Species == "scarlet"|
           Parasite.Species == "brood"|
           Parasite.Species == "pasteuria") %>%
  filter(AUC.single > 0)

colnames(host_auc_prev)[2] <- "Host.Species"

# eigenvector centrality scores for host/parasites of interest
host_eign_prev <- lake_diff_edges %>%
  filter(Weight.Method == "prevalence" &
           Metric == "eigenvector" &
           Parasite.Species != "metsch") 
  

testdf <- full_join(host_eign_prev, host_auc_prev)

testdf <- testdf %>% na.omit()

testdf$Parasite.Species <- as.factor(testdf$Parasite.Species)
testdf$Host.Species <- as.factor(testdf$Host.Species)

p1 <- ggplot(data = testdf,
             aes(x = AUC.single, y = Centrality, color = Host.Species)) +
  geom_point(size = 3, alpha = 0.5) +
  #geom_smooth(se = FALSE) +
  facet_grid(Host.Species~Parasite.Species, scales = "free")
p1


p2 <- ggplot(data = testdf,
             aes(x = log(AUC.single + 1), y = Centrality, color = Host.Species)) +
  geom_point(size = 3, alpha = 0.5) +
  #geom_smooth(se = FALSE) +
  facet_grid(Host.Species~Parasite.Species, scales = "free")
p2



# host species single AUC
host_auc_prev_lake <- lakeyear_auc_shared_prev %>%
  filter(Focal.Host == Other.Host) %>%
  group_by(Lake, Focal.Host, Parasite.Species) %>%
  summarise(AUC.single = sum(AUC.prev)) %>%
  filter(AUC.single > 0)





# (E) Prevalence epidemic size for brood, scarlet, and pasteuria in each, individual host species
epidemic_size <- lakeyear_auc_shared_prev %>%
  filter(Focal.Host == Other.Host) %>%
  filter(Focal.Host != "mendotae") %>%
  filter(AUC.prev > 0) %>%
  group_by(Lake, Parasite.Species, Focal.Host) %>% 
  summarise(AUC.avg = sum(AUC.prev)) %>% # NEW!!! summing across lakes
  mutate(AUC.avg.log = log(AUC.avg + 1)) %>%
  select(Lake, Parasite.Species, Focal.Host, AUC.avg, AUC.avg.log)
epidemic_size <- as.data.frame(epidemic_size)

epi_size_plot <- ggplot() +
  geom_boxplot(data = epidemic_size, aes(x = Focal.Host, y = AUC.avg), outlier.shape = NA,
               fill = "#f0f0f0") + 
  geom_jitter(data = epidemic_size, aes(x = Focal.Host,y = AUC.avg),
              width = .1, height = 0, size = 2, alpha = 0.4) +
  ylab("Single host epidemic size  (AUC of prevalence)") +
  xlab("") +
  ggtitle("") +
  facet_wrap(~Parasite.Species, scales = "free") + 
  coord_flip() 
epi_size_plot

epi_main_paras <- ggplot() +
  geom_boxplot(data = (epidemic_size %>% filter(Parasite.Species == "brood"|
                                                  Parasite.Species == "scarlet"|
                                                  Parasite.Species == "pasteuria")),
               aes(x = Focal.Host, y = AUC.avg), outlier.shape = NA,
               fill = "#f0f0f0") + 
  geom_jitter(data = (epidemic_size %>% filter(Parasite.Species == "brood"|
                                                 Parasite.Species == "scarlet"|
                                                 Parasite.Species == "pasteuria")), 
              aes(x = Focal.Host, y = AUC.avg),
              width = .1, height = 0, size = 2, alpha = 0.4) +
  ylab("Single host epidemic size (AUC of prevalence)") +
  xlab("") +
  ggtitle("") +
  facet_wrap(~Parasite.Species, scales = "free") + 
  coord_flip() 
epi_main_paras


kruskal.test(AUC.avg ~ Focal.Host, data = epidemic_size %>% filter(Parasite.Species == "scarlet"))

# stats for single epidemic AUC of prevalence
scar_epi_inf_kw <- kruskal.test(AUC.log ~ Focal.Host, data = scarlet_prev_epi_inf)
brood_epi_inf_kw <- kruskal.test(AUC.log ~ Focal.Host, data = brood_prev_epi_inf)
past_epi_inf_kw <- kruskal.test(AUC.log ~ Focal.Host, data = past_prev_epi_inf)

brood_con_i <- conover.test(brood_prev_epi_inf$AUC.log, brood_prev_epi_inf$Focal.Host, kw = TRUE, 
                            method = "bonferroni", altp = TRUE)
past_con_i <- conover.test(past_prev_epi_inf$AUC.log, past_prev_epi_inf$Focal.Host, kw = TRUE, 
                           method = "bonferroni", altp = TRUE)

# (F) Infecteds epidemic size for brood, scarlet, and pasteuria in each, individual host species
epidemic_size_inf <- lakeyear_auc_shared_prevdens %>%
  filter(Focal.Host == Other.Host) %>%
  filter(Focal.Host != "mendotae") %>%
  filter(AUC.prev > 0) %>%
  group_by(Parasite.Species, Focal.Host) %>% 
  mutate(AUC.log = log(AUC.prev + 1)) %>%
  select(Parasite.Species, Focal.Host, AUC.prev, AUC.log)
epidemic_size_inf <- as.data.frame(epidemic_size_inf)

epi_inf_plot <- ggplot() +
  geom_boxplot(data = epidemic_size_inf, aes(x = Focal.Host, y = AUC.log), outlier.shape = NA,
               fill = "#f0f0f0") + 
  geom_jitter(data = epidemic_size_inf, aes(x = Focal.Host,y = AUC.log),
              width = .1, height = 0, size = 2, alpha = 0.4) +
  ylab("Single host epidemic size  (log AUC of infected host density)") +
  xlab("") +
  ggtitle("") +
  facet_wrap(~Parasite.Species, scales = "free") + 
  coord_flip() 
epi_inf_plot

