# centrality-lake-figs

rm(list=ls())

setwd("~/Dropbox (University of Michigan)/Manuscripts/McIntireetal_Networks/Data&Code/daphnia-parasite-networks/2023/networks")
setwd("C:/Users/krism/Desktop/networks/networks")

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
lakeyear_auc_density <- read.csv("lakeyear_auc_density.csv", header = TRUE)
lakeyear_auc_shared_prevdens_geo <- read.csv("lakeyear_auc_shared_prevdens_geo.csv", header = TRUE)
lakeyear_auc_shared_total_density <- read.csv("lakeyear_auc_shared_total_density.csv", header = TRUE)

lake_auc_density <- lakeyear_auc_density %>%
  group_by(Lake, Host.Species) %>%
  summarise(AUC.density = sum(AUC.density)) %>%
  mutate(Log.AUC.density = log(AUC.density + 1))

# EPIDEMIC OVERLAP--INFECTED HOST DENSITY (geometric) ===========================================================

#arrange data for make_mat_list (don't log-transform version)
auc_shared <- lakeyear_auc_shared_prevdens_geo %>%
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
lakes <- length(lake_vec)
paras <- length(parasite_vec)
paralake <- lakes*paras

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
infecteds_geom_rank <- rank_lake_central

# COMPARING ACROSS METRICS ===========================================================
# make a new factor column for distinguishing data used to create edge weights
infecteds_geom_data <- infecteds_geom_rank %>% 
  mutate(Weight.Method = factor("infecteds_geom",
                                levels = c("infecteds_geom")))
common_parasite_vec <- c("brood", "scarlet", "pasteuria")


# PLOT: host centrality ================================================================

# (A) Eigenvector centrality--Pasteuria infecteds
xtext_vec <- seq(1, 5, by = 1)
mtext_vec <- c("b", "ab", "ab", "a", "a")
ytext_vec <- rep(1.1, 5)

subset_eigen_past <- infecteds_geom_data %>%
  filter(Metric == "eigenvector") %>%
  filter(Parasite.Species == "pasteuria")

eigen_past_inf_plot <- ggplot() +
  geom_boxplot(data = subset_eigen_past, aes(x = reorder(Host.Species, Centrality, FUN = median), y = Centrality),
               fill = "#377eb8", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(data = subset_eigen_past, aes(x = reorder(Host.Species, Centrality, FUN = median), 
                                       y = Centrality),
              width = .1, height = 0, size = 2, alpha = 0.4) +
  xlab("Host species") +
  ylab("Eigenvector centrality\n") +
  ggtitle(expression(italic("Infected host density"))) +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.0),
                     labels = c("0", "0.25", "0.5", "0.75", "1"),
                     limits = c(0, 1.125)) +
  ggtitle("P. ramosa") +
  annotate("text", x = xtext_vec, y = ytext_vec, label = c(mtext_vec)) +
  coord_flip() +
  theme_bw()
eigen_past_inf_plot

# (B) Eigenvector centrality--Brood infecteds
xtext_vec <- c(0.75, 1.75, 2.75, 3.75, 5, 6)
mtext_vec <- c("b", "b", "ab", "ab", "ab", "a")
ytext_vec <- rep(1.1, 6)

subset_eigen_brood <- infecteds_geom_data %>%
  filter(Metric == "eigenvector") %>%
  filter(Parasite.Species == "brood")

eigen_brood_inf_plot <- ggplot() +
  geom_boxplot(data = subset_eigen_brood, aes(x = reorder(Host.Species, Centrality, FUN = median), y = Centrality),
               fill = "#4daf4a", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(data = subset_eigen_brood, aes(x = reorder(Host.Species, Centrality, FUN = median), 
                                       y = Centrality),
              width = .1, height = 0, size = 2, alpha = 0.4) +
  xlab("") +
  ylab("Eigenvector centrality\n") +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.0),
                     labels = c("0", "0.25", "0.5", "0.75", "1"),
                     limits = c(0, 1.125)) +
  ggtitle("B.paedophthorum") +
  annotate("text", x = xtext_vec, y = ytext_vec, label = c(mtext_vec)) +
  coord_flip() +
  theme_bw()
eigen_brood_inf_plot

# (C) Eigenvector centrality--Scarlet infecteds
subset_eigen_spiro <- infecteds_geom_data %>%
  filter(Metric == "eigenvector") %>%
  filter(Parasite.Species == "scarlet")

eigen_scar_inf_plot <- ggplot() +
  geom_boxplot(data = subset_eigen_spiro, aes(x = reorder(Host.Species, Centrality, FUN = median), y = Centrality),
               fill = "#e41a1c", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(data = subset_eigen_spiro, aes(x = reorder(Host.Species, Centrality, FUN = median), 
                                       y = Centrality),
              width = .1, height = 0, size = 2, alpha = 0.4) +
  ylab("Eigenvector centrality\n") +
  xlab("") +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.0),
                     labels = c("0", "0.25", "0.5", "0.75", "1"),
                     limits = c(0, 1.125)) +
  ggtitle("S. cienkowskii") +
  coord_flip() + 
  theme_bw()
eigen_scar_inf_plot

centr_plot <- plot_grid(eigen_past_inf_plot,eigen_brood_inf_plot,
                        eigen_scar_inf_plot,
                        ncol = 3, labels = "AUTO")

pdf(file = "fig_centrality_geom_infdensity.pdf", width = 10, height = 3.5, useDingbats=F)
centr_plot
dev.off()


# STATS: centrality ================================================================
str(subset_eigen_past)
str(subset_eigen_brood)
str(subset_eigen_spiro)

kw_past <- kruskal.test(Centrality ~ Host.Species, data = subset_eigen_past)
kw_past
kw_brood <- kruskal.test(Centrality ~ Host.Species, data = subset_eigen_brood)
kw_brood
kw_spiro <- kruskal.test(Centrality ~ Host.Species, data = subset_eigen_spiro)
kw_spiro

comi_past <- conover.test(subset_eigen_past$Centrality, subset_eigen_past$Host.Species, kw = TRUE, 
                          method = "bonferroni", altp = TRUE)
comi_brood <- conover.test(subset_eigen_brood$Centrality, subset_eigen_brood$Host.Species, kw = TRUE, 
                          method = "bonferroni", altp = TRUE)



kw_result <- kruskal.test(Centrality ~ Host.Species, data = testdata)

coni_result <- conover.test(testdata$Centrality, testdata$Host.Species, kw = TRUE, 
                            method = "bonferroni", altp = TRUE)

chi_sq <- unname(kw_result$statistic) # get chi squared value from KW test
kw_pvalue <- kw_result$p.value # get overall p value from KW test

pwise_pvalue <- coni_result$altP.adjusted
pwise_comparison <- coni_result$comparisons





past_epi_inf_kw <- kruskal.test(AUC.log ~ Focal.Host, data = past_prev_epi_inf)

brood_con_i <- conover.test(brood_prev_epi_inf$AUC.log, brood_prev_epi_inf$Focal.Host, kw = TRUE, 
                            method = "bonferroni", altp = TRUE)







# Kruskal-Wallis test & Conover-Iman post hoc
edge_method_vec <- c("infecteds_geom")

raw_analyses_hosts <- data.frame(matrix(rep(NA, 120*6), nrow = 120)) 
names(raw_analyses_hosts) <- c("Parasite.Species", "Weight.Method", "KW.p", "KW.chisq", 
                               "Comparison", "Conover.p")
centrality_stats_conover <- compare_hosts_con(infecteds_geom_data, common_parasite_vec, edge_method_vec, raw_analyses_hosts)

centrality_kw_sig <- centrality_stats_conover %>%
  filter(KW.p < 0.05)

centrality_conover_sig <- centrality_stats_conover %>%
  filter(Conover.p < 0.05)

write.csv(centrality_conover_sig, "centrality_conover_sig.csv")


# PLOTS: AUC single epidemic infected host density ============================================

epidemic_size_inf <- lakeyear_auc_shared_prevdens_geo %>%
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


xtext_vec <- seq(1, 6, by = 1)
mtext_vec <- c("b", "ab", "ab", "ab", "a", "a")
ytext_vec <- rep(16.5, 6)

epi_past_inf_plot <- ggplot() +
  geom_boxplot(data = past_prev_epi_inf, aes(x = reorder(Focal.Host, AUC.log, FUN = median), y = AUC.log),
               fill = "#377eb8", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(data = past_prev_epi_inf, aes(x = reorder(Focal.Host, AUC.log, FUN = median), 
                                            y = AUC.log),
              width = .1, height = 0, size = 2, alpha = 0.4) +
  ylab("Epidemic size (log AUC) \ninfected host density") +
  xlab("") +
  ggtitle("P.ramosa") +
  scale_y_continuous(breaks = c(6, 8, 10, 12, 14, 16),
                     labels = c("6", "8", "10", "12", "14", "16"),
                     limits = c(6, 16.5)) +  
  annotate("text", x = xtext_vec, y = ytext_vec, label = c(mtext_vec)) +
  coord_flip() +
  theme_bw()
epi_past_inf_plot

xtext_vec <- seq(1, 6, by = 1)
mtext_vec <- c("c", "bc", "abc", "ab", "a", "abc")
ytext_vec <- rep(16, 6)

epi_brood_inf_plot <- ggplot() +
  geom_boxplot(data = brood_prev_epi_inf, aes(x = reorder(Focal.Host, AUC.log, FUN = median), y = AUC.log),
               fill = "#4daf4a", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(data = brood_prev_epi_inf, aes(x = reorder(Focal.Host, AUC.log, FUN = median), 
                                             y = AUC.log),
              width = .1, height = 0, size = 2, alpha = 0.4) +
  ylab("Epidemic size (log AUC) \ninfected host density") +
  xlab("") +
  ggtitle("B.paedophthorum") +
  scale_y_continuous(breaks = c(6, 8, 10, 12, 14, 16),
                     labels = c("6", "8", "10", "12", "14", "16"),
                     limits = c(6, 16.5)) +
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
  ylab("Epidemic size (log AUC) \ninfected host density") +
  xlab("") +
  ggtitle("S.cienkowskii") +
  scale_y_continuous(breaks = c(6, 8, 10, 12, 14, 16),
                     labels = c("6", "8", "10", "12", "14", "16"),
                     limits = c(6, 16.5)) +
  coord_flip() +
  theme_bw()
epi_scar_inf_plot


auc_plot <- plot_grid(epi_past_inf_plot,epi_brood_inf_plot, epi_scar_inf_plot,
                      ncol = 3, labels = "AUTO")

pdf(file = "fig_AUC_infhostdensity.pdf", width = 10, height = 3.5, useDingbats=F)
auc_plot
dev.off()

# stats for single epidemic AUC of prevalence
scar_epi_inf_kw <- kruskal.test(AUC.log ~ Focal.Host, data = scarlet_prev_epi_inf)
brood_epi_inf_kw <- kruskal.test(AUC.log ~ Focal.Host, data = brood_prev_epi_inf)
past_epi_inf_kw <- kruskal.test(AUC.log ~ Focal.Host, data = past_prev_epi_inf)


past_con_i <- conover.test(past_prev_epi_inf$AUC.log, past_prev_epi_inf$Focal.Host, kw = TRUE, 
                           method = "bonferroni", altp = TRUE)




# messing around space

