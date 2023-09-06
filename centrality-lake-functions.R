# lake centrality calculation functions


# "make_mat_list" creates a double list of host adjacency matrices at lake*parasite level
make_mat_list <- function(auc_shared, parasite_vec, lake_vec) {
  
  mat_list_list <- list()
  mat_list <- list()
  
  # filter and summarize auc data by parasite species
  for (i in 1:length(parasite_vec)) {
    paradat <- auc_shared %>%
      filter(Parasite.Species == parasite_vec[i]) 
    
    for (j in 1:length(lake_vec)) {
      lakedat <- paradat %>%
        filter(Lake == lake_vec[j]) %>%
        select(Focal.Host, Other.Host, Lake, Edge.Weight)
      
      dat <- as.data.frame(lakedat) # make into dataframe
      
      # turn data frame into matrix (using bipartite package function)
      dat_matrix <- frame2webs(dat, 
                               varnames = c("Focal.Host", "Other.Host", "Lake", "Edge.Weight"), 
                               type.out = "list", 
                               emptylist = FALSE) # if TRUE, drops empty rows/cols
      
      mat <- dat_matrix[[lake_vec[j]]]
      mat_list[[j]] <- mat # add matrix to list
      
    }
    names(mat_list) <- lake_vec
    mat_list_list[[i]] <- mat_list
    
  }
  return(mat_list_list) 
}

# "calc_graph_density" calculates graph density, node & edge number 
calc_graph_density <- function(auc_shared_mat, raw_graphdensity_data, parasite_vec, lake_vec, lake_auc_density) {
  
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
      
      # calculate metrics of interest
      num_vertex <- igraph::vcount(graph_trimmed)
      num_edge <- igraph::ecount(graph_trimmed)
      g_density <- igraph::edge_density(graph_trimmed, loops = TRUE)
      
      graph_vals <- c(lake_vec[j],
                      num_vertex, 
                      num_edge,
                      g_density)
      
      raw_graphdensity_data[x, 2:5] <- graph_vals
      x <- x + 1
      
    }
    par_vals <- c(parasite_vec[i])
    raw_graphdensity_data[y, 1] <- par_vals
    
    lakes <- length(lake_vec)
    y <- y + lakes
  }
  raw_graphdensity_data
}    

# "plot_graph_density" plots graph density data (w/ and w/o zeroes in the data)
plot_graph_density <- function(graph_density_data) {
  # reclass variables
  graph_density_data$Parasite.Species <- as.factor(graph_density_data$Parasite.Species)
  graph_density_data$Lake <- as.factor(graph_density_data$Lake)
  graph_density_data$Vertices <- as.integer(graph_density_data$Vertices)
  graph_density_data$Edges <- as.integer(graph_density_data$Edges)
  graph_density_data$Graph.Density <- as.numeric(graph_density_data$Graph.Density)
  
  # GRAPH DENSITY--boxplot with all data included
  graph_density <- ggplot(data = graph_density_data, 
                          aes(x = reorder(Parasite.Species, Graph.Density, FUN = median), 
                              y = Graph.Density)) +
    geom_boxplot(fill = "gray", outlier.shape = NA) +
    geom_jitter(data = graph_density_data, aes(x = reorder(Parasite.Species, Graph.Density, FUN = median),
                                               y = Graph.Density),
                width = .1, height = 0, size = 2, alpha = 0.4) +
    xlab("") +
    ylab("Graph density \n (total edges / total possible edges)") +
    ylim(0, 1) +
    ggtitle("Graph density (lake*parasite level networks)") +
     
    coord_flip()

  
  
  # GRAPH DENSITY--boxplot where we drop graphs with no edges 
  nonzerodata <- graph_density_data %>% filter(Graph.Density > 0) # filter data
  
  nonzero_graph_density <- ggplot(data = nonzerodata, 
                                  aes(x = reorder(Parasite.Species, Graph.Density, FUN = median), 
                                      y = Graph.Density)) +
    geom_boxplot(fill = "gray", outlier.shape = NA) +
    geom_jitter(data = nonzerodata, aes(x = reorder(Parasite.Species, Graph.Density, FUN = median),
                                        y = Graph.Density),
                width = .1, height = 0, size = 2, alpha = 0.4) +
    xlab("") +
    ylab("Graph density \n (total edges / total possible edges)") +
    ylim(0, 1) +
    ggtitle("Non-zero graph density (lake*parasite level networks)") +
     
    coord_flip()
  
  graph_density_plots <- list(graph_density, nonzero_graph_density)
  graph_density_plots
}

# "calc_centrality_bylake" calculates all the centrality values for each node in each network
calc_centrality_bylake <- function(auc_shared_mat, raw_centrality_data, parasite_vec, lake_vec) {
  
  for (i in 1:length(parasite_vec)) {
    paralist <- auc_shared_mat[[parasite_vec[i]]] # subset a list for each parasite
    
    for (j in 1:length(lake_vec)) {
      
      graph_thing <- graph_from_adjacency_matrix(paralist[[lake_vec[j]]],
                                                 mode = c("undirected"),
                                                 weighted = TRUE,
                                                 diag = FALSE) # create graph from list
      
      graph_thing <- igraph::delete.vertices(graph_thing, # drop the vertices without any connections
                                             which(igraph::degree(graph_thing) < 1))
      
      graph_thing_unw <- graph_from_adjacency_matrix(sign(paralist[[lake_vec[j]]]),
                                                     mode = c("undirected"),
                                                     weighted = TRUE,
                                                     diag = FALSE) # create graph from list
      
      graph_thing_unw <- igraph::delete.vertices(graph_thing_unw, # drop the vertices without any connections
                                                 which(igraph::degree(graph_thing_unw) < 1)) 
      
      # use graph_thing_unw to calculate degree (unweighted)
      degr_vals <- igraph::degree(graph_thing_unw)
      
      # record number of vertices & edges in graph (so we can drop appropriate ones later)
      num_vertex <- igraph::vcount(graph_thing)
      num_edge<- igraph::ecount(graph_thing)
      
      # need to calculate eigen_centrality & strength first, then we can change the edge weights...
      eign_vals <- igraph::eigen_centrality(graph_thing, weights = E(graph_thing)$weight)
      strn_vals <- igraph::strength(graph_thing,
                                    loops = FALSE)
      
      E(graph_thing)$weight <- mean(E(graph_thing)$weight)/E(graph_thing)$weight
      
      # calculate centrality metrics
      btwn_vals <- igraph::betweenness(graph_thing)
      clos_vals <- igraph::closeness(graph_thing)
      
      # vectors for each centrality metric
      eign_vals <- c(eign_vals$vector, "Parasite.Species" = paste(parasite_vec[i]),
                     "Lake" = paste(lake_vec[j]), "Metric" = "eigenvector",
                     "Vertices" = paste(num_vertex), "Edges" = paste(num_edge))
      btwn_vals <- c(btwn_vals, "Parasite.Species" = paste(parasite_vec[i]),
                     "Lake" = paste(lake_vec[j]), "Metric" = "betweenness",
                     "Vertices" = paste(num_vertex), "Edges" = paste(num_edge))
      clos_vals <- c(clos_vals, "Parasite.Species" = paste(parasite_vec[i]),
                     "Lake" = paste(lake_vec[j]), "Metric" = "closeness",
                     "Vertices" = paste(num_vertex), "Edges" = paste(num_edge))
      degr_vals <- c(degr_vals, "Parasite.Species" = paste(parasite_vec[i]),
                     "Lake" = paste(lake_vec[j]), "Metric" = "degree",
                     "Vertices" = paste(num_vertex), "Edges" = paste(num_edge))
      strn_vals <- c(strn_vals, "Parasite.Species" = paste(parasite_vec[i]),
                     "Lake" = paste(lake_vec[j]), "Metric" = "strength",
                     "Vertices" = paste(num_vertex), "Edges" = paste(num_edge))
      
      # add vectors as rows to data frame
      raw_centrality_data <- bind_rows(raw_centrality_data, eign_vals, btwn_vals,
                                       clos_vals, degr_vals, strn_vals)
      
      if (j == 1 & i == 1) { # drop first row (it's blank to make the original data frame)
        raw_centrality_data <- raw_centrality_data[-1, ]
      }
    }
  }
  raw_centrality_data
}    


# "calc_diff_eigen" calculates weighted and unweighted eigenvector centrality scores
calc_diff_eigen <- function(auc_shared_mat, raw_centrality_data, parasite_vec, lake_vec) {
  
  for (i in 1:length(parasite_vec)) {
    paralist <- auc_shared_mat[[parasite_vec[i]]] # subset a list for each parasite
    
    for (j in 1:length(lake_vec)) {
      
      graph_thing <- graph_from_adjacency_matrix(paralist[[lake_vec[j]]],
                                                 mode = c("undirected"),
                                                 weighted = TRUE,
                                                 diag = FALSE) # create graph from list
      
      graph_thing <- igraph::delete.vertices(graph_thing, # drop the vertices without any connections
                                             which(igraph::degree(graph_thing) < 1))
      
      graph_thing_unw <- graph_from_adjacency_matrix(sign(paralist[[lake_vec[j]]]),
                                                     mode = c("undirected"),
                                                     weighted = TRUE,
                                                     diag = FALSE) # create graph from list
      
      graph_thing_unw <- igraph::delete.vertices(graph_thing_unw, # drop the vertices without any connections
                                                 which(igraph::degree(graph_thing_unw) < 1)) 
      
      
      # unweighted eigen
      eign_unw <- igraph::eigen_centrality(graph_thing_unw)
      
      # weighted eigen
      eign_wtd <- igraph::eigen_centrality(graph_thing, directed = FALSE, weights = E(graph_thing)$weight)
      
      # record number of vertices & edges in graph (so we can drop appropriate ones later)
      num_vertex <- igraph::vcount(graph_thing)
      num_edge<- igraph::ecount(graph_thing)
      
      # vectors for each centrality metric
      eign_unw_vals <- c(eign_unw$vector, "Parasite.Species" = paste(parasite_vec[i]),
                         "Lake" = paste(lake_vec[j]), "Metric" = "eigenvector_unw",
                         "Vertices" = paste(num_vertex), "Edges" = paste(num_edge))
      
      eign_wtd_vals <- c(eign_wtd$vector, "Parasite.Species" = paste(parasite_vec[i]),
                         "Lake" = paste(lake_vec[j]), "Metric" = "eigenvector_wtd",
                         "Vertices" = paste(num_vertex), "Edges" = paste(num_edge))
      
      # add vectors as rows to data frame
      raw_centrality_data <- bind_rows(raw_centrality_data, eign_unw_vals, eign_wtd_vals)
      
      if (j == 1 & i == 1) { # drop first row (it's blank to make the original data frame)
        raw_centrality_data <- raw_centrality_data[-1, ]
      }
    }
  }
  raw_centrality_data
}    


# "organize_centrality" cleans/organizes the raw centrality data
organize_centrality <- function(centrality_bylake) {
  
  # clean the data 
  lake_central <- centrality_bylake %>% 
    select(-ambigua, - mendotae) %>%
    filter(Vertices > 2) %>%
    gather(Host.Species, Centrality, ceriodaphnia:retrocurva) %>%
    na.omit()
  
  # change class of each column (all were converted to chr b/c of missing values)
  lake_central$Vertices <- as.integer(lake_central$Vertices) # integer columns
  lake_central$Edges <- as.integer(lake_central$Edges)
  lake_central$Lake <- as.factor(lake_central$Lake) # factor columns
  lake_central$Parasite.Species <- as.factor(lake_central$Parasite.Species)
  lake_central$Metric <- as.factor(lake_central$Metric)
  lake_central$Host.Species <- as.factor(lake_central$Host.Species)
  lake_central$Centrality <- as.numeric(lake_central$Centrality) # numeric columns
  
  # rank transform the centrality scores
  rank_lake_central <- lake_central %>%
    group_by(Parasite.Species, Lake, Metric) %>%
    mutate(Centrality.Rank = min_rank(-Centrality),
           Centrality.PercRank = percent_rank(-Centrality)) %>%
    group_by(Parasite.Species, Host.Species, Metric) %>%
    mutate(cent_mean = mean(Centrality))
  rank_lake_central <- as.data.frame(rank_lake_central) # make sure it's a df
  
  rank_lake_central # return this 
}


# "plot_raw_centrality" creates boxplots of the raw data, separately for each metric
plot_raw_centrality <- function(rank_lake_central, metric_vec, edge_method) {
  j <- 1 # dummy variable for list of plots  
  
  #raw_plots <- list()
  raw_plots <- vector("list", length(metric_vec))
  names(raw_plots) <- metric_vec
  
  for (i in 1:length(metric_vec)) {
    subset_metric <- rank_lake_central %>%
      filter(Metric == metric_vec[i])
    
    # plot the raw data
    raw_plot <- ggplot() +
      # geom_violin(data = subset_metric, aes(x = reorder(Host.Species, Centrality, FUN = mean), y = Centrality)) +
      # geom_jitter(data = subset_metric, aes(x = reorder(Host.Species, Centrality, FUN = mean), y = Centrality), 
      #             width = .1, height = 0, alpha = 0.4, size = 2) +
      # geom_point(data = subset_metric, aes(x = reorder(Host.Species, Centrality, FUN = mean), y = cent_mean), 
      #            shape = 18, color = "#1f78b4", size = 3) +
      geom_boxplot(data = subset_metric, aes(x = reorder(Host.Species, Centrality, FUN = mean), y = Centrality),
                   fill = "gray", outlier.shape = NA) +
      geom_jitter(data = subset_metric, aes(x = reorder(Host.Species, Centrality, FUN = mean), 
                                            y = Centrality),
                  width = .1, height = 0, size = 2, alpha = 0.4) +
      ylab(paste(metric_vec[i], "Centrality (raw scores)")) +
      xlab("") +
      ggtitle(paste(metric_vec[i], "(lake*parasite level)", "\n Edge =", edge_method)) +
      facet_wrap(~Parasite.Species, ncol = 2) + 
       
      coord_flip() 
    
    raw_plots[[j]] <- raw_plot
    j <- j + 1
  }
  raw_plots
}

# "plot_raw_centrality" creates boxplots of the raw data, separately for each metric
plot_raw_together <- function(rank_lake_central, edge_method) {
  raw_plot <- ggplot() +
    geom_boxplot(data = rank_lake_central, aes(x = reorder(Host.Species, Centrality, FUN = mean), y = Centrality),
                 fill = "gray", outlier.shape = NA) +
    geom_jitter(data = rank_lake_central, aes(x = reorder(Host.Species, Centrality, FUN = mean), 
                                              y = Centrality),
                width = .1, height = 0, size = 2, alpha = 0.4) +
    ylab("Centrality (raw scores)") +
    xlab("") +
    ggtitle(paste("(lake*parasite level)", "Edge =", edge_method)) +
    facet_grid(Parasite.Species~Metric, scales = "free") + 
     
    coord_flip() 
  
  raw_plot
}

# "plot_rank_centrality" plots the rank-transformed centrality data 
plot_rank_dot <- function(rank_lake_central, edge_method) {
  
  # get means and medians for each centrality
  mean_centrality <- rank_lake_central %>%
    group_by(Host.Species, Parasite.Species, Metric) %>%
    mutate(meanperc = mean(Centrality.PercRank)) 
  mean_centrality <- as.data.frame(mean_centrality)
  
  # get tallies of the number of times of each centrality rank
  count_centrality <- mean_centrality %>%
    group_by(Parasite.Species, Host.Species, Metric, Centrality.PercRank) %>%
    mutate(count = n())
  count_centrality <- as.data.frame(count_centrality)
  
  #rank_plots <- list()
  rank_plots <- vector("list", length(metric_vec))
  names(rank_plots) <- metric_vec
  
  # change parasite order in appropriate dfs
  count_centrality$Parasite.Species_f <- factor(count_centrality$Parasite.Species, 
                                                levels = c("brood", "scarlet", "pasteuria", "metsch"))
  mean_centrality$Parasite.Species_f <- factor(mean_centrality$Parasite.Species, 
                                               levels = c("brood", "scarlet", "pasteuria", "metsch"))
  count_centrality$Metric_f <- factor(count_centrality$Metric, 
                                      levels = c("betweenness", "degree", "strength", "closeness", "eigenvector"))
  mean_centrality$Metric_f <- factor(mean_centrality$Metric, 
                                     levels = c("betweenness", "degree", "strength", "closeness", "eigenvector"))
  
  # rank transformed centrality data plot (size of circles represents # times host had that centrality rank)
  rank_plot <- ggplot() +
    geom_point(data = count_centrality, aes(x = reorder(Host.Species, Centrality.PercRank, FUN = mean), 
                                            y = Centrality.PercRank, size = count), shape = 1, color = "#1f78b4") +
    geom_point(data = mean_centrality, aes(x = Host.Species, y = meanperc), 
               shape = 18, color = "#1f78b4", size = 3) +
    scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) +
    scale_size_continuous(name = "Number of \nobservations",range = c(1, 6)) +
    ylab("Centrality (Proportional Rank)") +
    xlab("") +
    ggtitle(paste("(lake*parasite level) Edge =", edge_method)) +
    facet_grid(Parasite.Species_f~Metric_f) + 
    
    coord_flip() +
    theme(legend.position="bottom")
  
  rank_plot
}

# "plot_rank_centrality" plots the rank-transformed centrality data 
plot_rank_box <- function(rank_lake_central, edge_method) {
  rank_plot <- ggplot() +
    geom_boxplot(data = rank_lake_central, aes(x = reorder(Host.Species, Centrality.PercRank, FUN = mean), y = Centrality.PercRank),
                 fill = "gray", outlier.shape = NA) +
    geom_jitter(data = rank_lake_central, aes(x = reorder(Host.Species, Centrality.PercRank, FUN = mean), 
                                                  y = Centrality.PercRank),
                width = .1, height = 0, size = 2, alpha = 0.4) +
    ylab("Centrality (rank scores)") +
    xlab("") +
    ggtitle(paste("(lake*parasite level) Edge =", edge_method)) +
    facet_grid(Parasite.Species~Metric, scales = "free") + 
     
    coord_flip() 
  
  rank_plot
}

# "plot_edge_comparison" plots boxplots of eigenvector centrality across multiple edge weight options
plot_edge_comparison_raw <- function(lake_diff_edges, common_parasite_vec) {
  j <- 1 # dummy variable for list of plots  
  
  # eigenvector only for now
  eigen_lake_diff <- lake_diff_edges %>% filter(Metric == "eigenvector")
  
  edge_plots <- vector("list", length(common_parasite_vec))
  names(edge_plots) <- common_parasite_vec
  
  for (i in 1:length(common_parasite_vec)) {
    subset_parasite <- eigen_lake_diff %>%
      filter(Parasite.Species == common_parasite_vec[i])
    
    edge_plot <- ggplot() +
      geom_boxplot(data = subset_parasite, aes(x = reorder(Host.Species, Centrality, FUN = mean), y = Centrality),
                   fill = "gray", outlier.shape = NA) +
      geom_jitter(data = subset_parasite, aes(x = reorder(Host.Species, Centrality, FUN = mean), 
                                                        y = Centrality),
                  width = .1, height = 0, size = 2, alpha = 0.4) +
      ylab("Centrality (raw scores)") +
      xlab("") +
      ggtitle("Comparing diff. edge weights \nLake level networks") +
      facet_grid(Weight.Method~Parasite.Species) + 
       
      coord_flip() 
    edge_plot
    
    
    edge_plots[[j]] <- edge_plot
    j <- j + 1
  }
  edge_plots
}


# "plot_edge_comparison" plots boxplots of eigenvector centrality across multiple edge weight options
plot_edge_comparison_rank <- function(lake_diff_edges, common_parasite_vec) {
  j <- 1 # dummy variable for list of plots  
  
  # eigenvector only for now
  eigen_lake_diff <- lake_diff_edges %>% filter(Metric == "eigenvector")
  
  edge_plots <- vector("list", length(common_parasite_vec))
  names(edge_plots) <- common_parasite_vec
  
  for (i in 1:length(common_parasite_vec)) {
    subset_parasite <- eigen_lake_diff %>%
      filter(Parasite.Species == common_parasite_vec[i])
    
    edge_plot <- ggplot() +
      geom_boxplot(data = subset_parasite, aes(x = reorder(Host.Species, -Centrality.PercRank, FUN = mean), y = Centrality.PercRank),
                   fill = "gray", outlier.shape = NA) +
      geom_jitter(data = subset_parasite, aes(x = reorder(Host.Species, -Centrality.PercRank, FUN = mean), 
                                              y = Centrality.PercRank),
                  width = .1, height = 0, size = 2, alpha = 0.4) +
      ylab("Centrality (rank percent)") +
      xlab("") +
      ggtitle("Comparing diff. edge weights \nLake level networks") +
      facet_grid(Weight.Method~Parasite.Species) + 
       
      coord_flip() 
    edge_plot
    
    
    edge_plots[[j]] <- edge_plot
    j <- j + 1
  }
  edge_plots
}

# "compare_hosts" calculates kruskal-wallis and dunn's post hoc test for centrality data
compare_hosts <- function(lake_diff_edges, common_parasite_vec, edge_method_vec, raw_analyses_hosts) {
  x <- 1

  for (i in 1:length(common_parasite_vec)) {
    subset_parasite <- lake_diff_edges %>%
      filter(Parasite.Species == common_parasite_vec[i] & Metric == "eigenvector")
    
    for (j in 1:length(edge_method_vec)) {
      testdata <- subset_parasite %>% filter(Weight.Method == edge_method_vec[j])
      
      kw_result <- kruskal.test(Centrality ~ Host.Species, data = testdata)
      dunn_result <- dunn.test(testdata$Centrality, testdata$Host.Species, kw = TRUE, 
                               method = "bonferroni", altp = TRUE)
      
      chi_sq <- unname(kw_result$statistic) # get chi squared value from KW test
      kw_pvalue <- kw_result$p.value # get overall p value from KW test
      
      pwise_pvalue <- dunn_result$altP.adjusted
      pwise_comparison <- dunn_result$comparisons
      
      y <- length(pwise_comparison)
      z <- x + y - 1
      
      kw_pvalue_rep <- rep(kw_pvalue, y)
      chi_sq_rep <- rep(chi_sq, y)
      
      raw_analyses_hosts[x:z, 1] <- common_parasite_vec[i]
      raw_analyses_hosts[x:z, 2] <- edge_method_vec[j]
      raw_analyses_hosts[x:z, 3] <- kw_pvalue_rep
      raw_analyses_hosts[x:z, 4] <- chi_sq_rep
      raw_analyses_hosts[x:z, 5] <- pwise_comparison
      raw_analyses_hosts[x:z, 6] <- pwise_pvalue
        
      x <- x + y
    }
  }
  raw_analyses_hosts
}


# "compare_parasites" calculates kruskal-wallis and dunn's post hoc test for graph density data
compare_parasites <- function(nonzero_graph_density, raw_analyses_parasites) {
  
  kw_result <- kruskal.test(Graph.Density ~ Parasite.Species, data = nonzero_graph_density)
  dunn_result <- dunn.test(nonzero_graph_density$Graph.Density, nonzero_graph_density$Parasite.Species,
                           kw = TRUE, method = "bonferroni", altp = TRUE)
  
  chi_sq <- unname(kw_result$statistic) # get chi squared value from KW test
  kw_pvalue <- kw_result$p.value # get overall p value from KW test
  
  pwise_pvalue <- dunn_result$altP.adjusted
  pwise_comparison <- dunn_result$comparisons
  
  y <- length(pwise_comparison)
  kw_pvalue_rep <- rep(kw_pvalue, y)
  chi_sq_rep <- rep(chi_sq, y)
  
  raw_analyses_parasites[, 1] <- kw_pvalue_rep
  raw_analyses_parasites[, 2] <- chi_sq_rep
  raw_analyses_parasites[, 3] <- pwise_comparison
  raw_analyses_parasites[, 4] <- pwise_pvalue

  raw_analyses_parasites
}



# "compare_hosts_con" calculates kruskal-wallis and conover-iman post hoc test for centrality data
compare_hosts_con <- function(lake_diff_edges, common_parasite_vec, edge_method_vec, raw_analyses_hosts) {
  x <- 1
  
  for (i in 1:length(common_parasite_vec)) {
    subset_parasite <- lake_diff_edges %>%
      filter(Parasite.Species == common_parasite_vec[i] & Metric == "eigenvector")
    
    for (j in 1:length(edge_method_vec)) {
      testdata <- subset_parasite %>% filter(Weight.Method == edge_method_vec[j])
      
      kw_result <- kruskal.test(Centrality ~ Host.Species, data = testdata)

      coni_result <- conover.test(testdata$Centrality, testdata$Host.Species, kw = TRUE, 
                               method = "bonferroni", altp = TRUE)
      
      chi_sq <- unname(kw_result$statistic) # get chi squared value from KW test
      kw_pvalue <- kw_result$p.value # get overall p value from KW test
      
      pwise_pvalue <- coni_result$altP.adjusted
      pwise_comparison <- coni_result$comparisons
      
      y <- length(pwise_comparison)
      z <- x + y - 1
      
      kw_pvalue_rep <- rep(kw_pvalue, y)
      chi_sq_rep <- rep(chi_sq, y)
      
      raw_analyses_hosts[x:z, 1] <- common_parasite_vec[i]
      raw_analyses_hosts[x:z, 2] <- edge_method_vec[j]
      raw_analyses_hosts[x:z, 3] <- kw_pvalue_rep
      raw_analyses_hosts[x:z, 4] <- chi_sq_rep
      raw_analyses_hosts[x:z, 5] <- pwise_comparison
      raw_analyses_hosts[x:z, 6] <- pwise_pvalue
      
      x <- x + y
    }
  }
  raw_analyses_hosts
}


# new version for all the metrics?
# plot_edge_comparison_rank <- function(lake_diff_edges, common_parasite_vec, metric_vec) {
#   j <- 1 # dummy variable for list of plots  
#   
#   # eigenvector only for now
#   eigen_lake_diff <- lake_diff_edges %>% filter(Metric == "eigenvector")
#   
#   edge_plot_list <- vector("list", length(common_parasite_vec))
#   names(edge_plots) <- common_parasite_vec
#   
#   edge_plot_list_list <- vector("list", length(metric_vec))
#   names(edge_plots) <- metric_vec
#   
#   for (i in 1:length(common_parasite_vec)) {
#     subset_parasite <- eigen_lake_diff %>%
#       filter(Parasite.Species == common_parasite_vec[i])
#     
#     edge_plot <- ggplot() +
#       geom_boxplot(data = subset_parasite, aes(x = reorder(Host.Species, -Centrality.PercRank, FUN = mean), y = Centrality.PercRank),
#                    fill = "gray", outlier.shape = NA) +
#       geom_jitter(data = subset_parasite, aes(x = reorder(Host.Species, -Centrality.PercRank, FUN = mean), 
#                                               y = Centrality.PercRank),
#                   width = .1, height = 0, size = 2, alpha = 0.4) +
#       ylab("Centrality (rank percent)") +
#       xlab("") +
#       ggtitle("Comparing diff. edge weights \nLake level networks") +
#       facet_grid(Weight.Method~Parasite.Species) + 
#        
#       coord_flip() 
#     edge_plot
#     
#     
#     edge_plots[[j]] <- edge_plot
#     j <- j + 1
#   }
#   edge_plots
# }
# 


