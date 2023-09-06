# breaking down graph density, host breadth, and epidemic overlap

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

# new function to break things down more

# "calc_graph_density" calculates graph density, node & edge number 
calc_graph_density_2 <- function(auc_shared_mat, raw_graphdensity_data, parasite_vec, lake_vec, lake_auc_density) {
  
  # set up dummy variables to control position of data entry
  x <- 1
  y <- 1:15
  
  
  for (i in 1:length(parasite_vec)) {
    
    paralist <- auc_shared_mat[[parasite_vec[i]]]
    
    for (j in 1:length(lake_vec)) {
      lakeparamat <- paralist[[lake_vec[j]]]
      
      graph_thing <- graph_from_adjacency_matrix(lakeparamat,
                                                 mode = c("undirected"),
                                                 weighted = TRUE,
                                                 diag = TRUE)
      #NOTE: diag = TRUE here b/c we care about single host infections too
      
      densities <- lake_auc_density %>% # filter density data for lake*year
        filter(Lake == lake_vec[j])
      
      # get vector of host names
      absent_hosts <- densities$Host.Species[densities$AUC.density == 0]
      
      # drop absent hosts from graph
      graph_trimmed <- delete_vertices(graph_thing, paste(absent_hosts))
      
      # drop loops from graph
      graph_unloop <- simplify(graph_trimmed)
      
      # calculate metrics of interest
      num_vertex <- igraph::vcount(graph_trimmed) # total nodes
      num_edge <- igraph::ecount(graph_trimmed) # total realized edges, including loops
      g_density <- igraph::edge_density(graph_trimmed, loops = TRUE) # igraph graph density, with loops
      num_nonloop <- igraph::ecount(graph_unloop) # number realized edges, without loops
      nonloop_density <- igraph::edge_density(graph_unloop, loops = FALSE) # igraph graph density, without loops
      complete_noloop <- num_vertex*(num_vertex - 1)/2 
      
      graph_vals <- c(lake_vec[j],
                      num_vertex, 
                      num_edge,
                      num_nonloop,
                      complete_noloop,
                      g_density,
                      nonloop_density)
      
      raw_graphdensity_data[x, 2:8] <- graph_vals
      x <- x + 1
      
    }
    par_vals <- c(parasite_vec[i])
    raw_graphdensity_data[y, 1] <- par_vals
    
    lakes <- length(lake_vec)
    y <- y + lakes
  }
  raw_graphdensity_data
}    


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



# density of each host species in each lake
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

parasite_vec <- levels(auc_shared$Parasite.Species) # set up vectors of names for loop
host_vec <- levels(auc_shared$Focal.Host)
lake_vec <- levels(auc_shared$Lake)
lakes <- length(lake_vec)
paras <- length(parasite_vec)
paralake <- lakes*paras

# make list of host adjacency matrices
auc_shared_mat <- make_mat_list(auc_shared, parasite_vec, lake_vec)  
names(auc_shared_mat) <- parasite_vec

# calculate and plot graph density
raw_graphdensity_data <- data.frame(matrix(rep(NA, paralake*8), nrow = paralake))
names(raw_graphdensity_data) <- c("Parasite.Species", "Lake","Vertices", "Edges", 
                                  "NonLoops", "Complete","Graph.Density", "Overlap.Density")

graph_density_data <- calc_graph_density_2(auc_shared_mat, raw_graphdensity_data, parasite_vec, lake_vec, lake_auc_density)

# graph_density_data$Parasite.Species <- as.factor(graph_density_data$Parasite.Species)
# graph_density_data$Lake <- as.factor(graph_density_data$Lake)
# graph_density_data$Vertices <- as.integer(graph_density_data$Vertices)
# graph_density_data$Edges <- as.integer(graph_density_data$Edges)
# graph_density_data$NonLoops <- as.integer(graph_density_data$NonLoops)
# graph_density_data$Complete <- as.integer(graph_density_data$Complete)
# graph_density_data$Graph.Density <- as.numeric(graph_density_data$Graph.Density)
# graph_density_data$Overlap.Density <- as.numeric(graph_density_data$Overlap.Density)
# 
# str(test)
# 
# test <- graph_density_data %>%
#   filter(Parasite.Species == "brood" & Lake == "Appleton") %>%
#   mutate(Loops = Edges - NonLoops) %>%
#   mutate(HostBreadth = Loops/Vertices) %>%
#   mutate(fullgraphdensityA = (2*NonLoops + Loops)/(Vertices*(Vertices - 1) + Vertices)) %>%
#   mutate(fullgraphdensityB = (2*NonLoops + 2*Loops)/(Vertices*(Vertices - 1) + 2*Vertices)) %>%
#   mutate(nonloopgradphdensity = 2*NonLoops/(Vertices*(Vertices - 1)))


graph_density_plots <- plot_graph_density(graph_density_data)
names(graph_density_plots) <- c("all_data", "nonzero_data")


# GRAPH DENSITY ======================================

# calculate graph density stats
nonzero_graph_density <- graph_density_data %>% filter(Graph.Density > 0) # filter data
nonzero_graph_density$Parasite.Species <- as.factor(nonzero_graph_density$Parasite.Species)
nonzero_graph_density$Lake <- as.factor(nonzero_graph_density$Lake)
nonzero_graph_density$Vertices <- as.integer(nonzero_graph_density$Vertices)
nonzero_graph_density$Edges <- as.integer(nonzero_graph_density$Edges)
nonzero_graph_density$NonLoops <- as.integer(nonzero_graph_density$NonLoops)
nonzero_graph_density$Complete <- as.integer(nonzero_graph_density$Complete)
nonzero_graph_density$Graph.Density <- as.numeric(nonzero_graph_density$Graph.Density)
nonzero_graph_density$Overlap.Density <- as.numeric(nonzero_graph_density$Overlap.Density)

# add number of epidemic overlaps as a measure
str(nonzero_graph_density)

gdensity_clean <- nonzero_graph_density %>%
  mutate(Loops = Edges - NonLoops) %>%
  mutate(HostBreadth = Loops/Vertices) %>%
  mutate(Graph.Density.new = (2*NonLoops + Loops)/(Vertices*(Vertices - 1) + Vertices))

str(gdensity_clean)

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

kw_result <- kruskal.test(Graph.Density.new ~ Parasite.Species, data = gdensity_clean)
dunn_result <- dunn.test(gdensity_clean$Graph.Density.new, gdensity_clean$Parasite.Species,
                         kw = TRUE, method = "bonferroni", altp = TRUE)
coniman_result <- conover.test(gdensity_clean$Graph.Density.new, gdensity_clean$Parasite.Species, kw = TRUE, 
                               method = "bonferroni", altp = TRUE)

# PLOT 1: graph density ================================================================

# Graph density with conover-iman post-hoc
max_gdens <- gdensity_clean %>%
  group_by(Parasite.Species) %>%
  summarise(max = max(Graph.Density.new)) %>%
  arrange(max)

xtext_vec <- seq(1, 7, by = 1)
mtext_vec <- c("bc", "c", "bc", "bc", "ab", "ab", "a")
ytext_vec <- c(max_gdens$max + 0.255)


gdensity_clean$Parasite.Species <- str_replace(gdensity_clean$Parasite.Species, "scarlet", "spiro")
gdensity_clean$Parasite.Species <- as.factor(gdensity_clean$Parasite.Species)


gdensity_plot <- ggplot(data = gdensity_clean, 
                        aes(x = reorder(Parasite.Species, Graph.Density, FUN = median), 
                            y = Graph.Density.new)) +
  geom_boxplot(fill = "#f0f0f0", outlier.shape = NA) +
  geom_jitter(data = gdensity_clean, aes(x = reorder(Parasite.Species, Graph.Density, FUN = median),
                                                y = Graph.Density.new),
              width = 0.1, height = 0, size = 2.5, alpha = 0.4) +
  xlab("") +
  scale_x_discrete(labels=c("G.vavrai", "L.obtusa","Spider", "M.bicuspidata","P.ramosa","S.cienkowskii","B.paedophthorum"))+
  ylab("Graph density") +
  ylim(0, 1) +
  ggtitle("") +
  annotate("text", x = xtext_vec, y = ytext_vec, label = c(mtext_vec)) +
  coord_flip() +
  theme_bw()
gdensity_plot

# host breadth with conover-iman post-hoc
max_gdens <- gdensity_clean %>%
  group_by(Parasite.Species) %>%
  summarise(max = max(HostBreadth)) %>%
  arrange(max)

xtext_vec <- seq(1, 7, by = 1)
mtext_vec <- c("bc", "c", "bc", "bc", "ab", "ab", "a")
ytext_vec <- c(max_gdens$max + 0.125)
ytext_vec[7] <- 0.9250000 # so it doesn't go off

hostbreadth_plot <- ggplot(data = gdensity_clean, 
                        aes(x = reorder(Parasite.Species, HostBreadth, FUN = median), 
                            y = HostBreadth)) +
  geom_boxplot(fill = "#f0f0f0", outlier.shape = NA) +
  geom_jitter(data = gdensity_clean, aes(x = reorder(Parasite.Species, HostBreadth, FUN = median),
                                         y = HostBreadth),
              width = 0.1, height = 0, size = 2.5, alpha = 0.4) +
  xlab("") +
  scale_x_discrete(labels=c("G.vavrai", "L.obtusa","Spider", "M.bicuspidata","P.ramosa","S.cienkowskii","B.paedophthorum"))+
  ylab("Host breadth") +
  ylim(0, 1) +
  ggtitle("") +
  annotate("text", x = xtext_vec, y = ytext_vec, label = c(mtext_vec)) +
  coord_flip() +
  theme_bw()
hostbreadth_plot

# epidemic overlaps with conover-iman post-hoc
max_gdens <- gdensity_clean %>%
  group_by(Parasite.Species) %>%
  summarise(max = max(Overlap.Density)) %>%
  arrange(max)

xtext_vec <- seq(1, 7, by = 1)
mtext_vec <- c("c", "c", "bc", "bc", "ab", "ab", "a")
ytext_vec <- c(max_gdens$max + 0.125)
ytext_vec <- ytext_vec[c(1:5, 7, 6)] # reorder

epioverlap_plot <- ggplot(data = gdensity_clean, 
                           aes(x = reorder(Parasite.Species, Overlap.Density, FUN = median), 
                               y = Overlap.Density)) +
  geom_boxplot(fill = "#f0f0f0", outlier.shape = NA) +
  geom_jitter(data = gdensity_clean, aes(x = reorder(Parasite.Species, Overlap.Density, FUN = median),
                                         y = Overlap.Density),
              width = 0.1, height = 0, size = 2.5, alpha = 0.4) +
  xlab("") +
  scale_x_discrete(labels=c("G.vavrai", "L.obtusa", "M.bicuspidata","Spider","P.ramosa","S.cienkowskii","B.paedophthorum"))+
  ylab("Cross species transmission") +
  ylim(0, 1) +
  ggtitle("") +
  annotate("text", x = xtext_vec, y = ytext_vec, label = c(mtext_vec)) +
  coord_flip() +
  theme_bw()
epioverlap_plot


para_plots <- plot_grid(hostbreadth_plot, epioverlap_plot, gdensity_plot,
                        ncol = 3, labels = "AUTO")

pdf(file = paste("fig_graphdensity_and_more.pdf"), width = 10, height = 4, useDingbats=F)
para_plots
dev.off()


str(gdensity_clean)



kw_result_A <- kruskal.test(Graph.Density ~ Parasite.Species, data = gdensity_clean)
kw_result_B <- kruskal.test(HostBreadth ~ Parasite.Species, data = gdensity_clean)
kw_result_C <- kruskal.test(Overlap.Density ~ Parasite.Species, data = gdensity_clean)

coniman_result_A <- conover.test(gdensity_clean$Graph.Density, gdensity_clean$Parasite.Species, kw = TRUE, 
                               method = "bonferroni", altp = TRUE)
coniman_result_B <- conover.test(gdensity_clean$HostBreadth, gdensity_clean$Parasite.Species, kw = TRUE, 
                                 method = "bonferroni", altp = TRUE)
coniman_result_C <- conover.test(gdensity_clean$Overlap.Density, gdensity_clean$Parasite.Species, kw = TRUE, 
                                 method = "bonferroni", altp = TRUE)





dunn_result <- dunn.test(nonzero_graph_density$Graph.Density, nonzero_graph_density$Parasite.Species,
                         kw = TRUE, method = "bonferroni", altp = TRUE)
coniman_result <- conover.test(nonzero_graph_density$Graph.Density, nonzero_graph_density$Parasite.Species, kw = TRUE, 
                               method = "bonferroni", altp = TRUE)



pdf(file = paste("fig_graphdensity.pdf"), width = 4.25, height = 4, useDingbats=F)
gdensity_plot
dev.off()

# get means (for references in paper)
mean_gdensity <- nonzero_graph_density %>% 
  group_by(Parasite.Species) %>%
  summarise(mn = mean(Graph.Density))



