require(igraph)

# Graph A ------------------------------------------------------------

# adjacency matrix values
adj_mat_A <- matrix(c(1,0,0,0,0,
                      0,1,0,0,0,
                      0,0,1,0,0,
                      0,0,0,1,0,
                      0,0,0,0,0),
                    nrow = 5,
                    byrow = TRUE) 

# add node names
dimnames(adj_mat_A) = list(c("A", "B", "C", "D", "E"), 
                           c("A", "B", "C", "D", "E")) 

p1 <- graph_from_adjacency_matrix(adj_mat_A,
                                  mode = "undirected",
                                  weighted = TRUE)

plot(p1, layout = layout.fruchterman.reingold,
     vertex.color = "#80b1d3",
     edge.color = "black",
     edge.width = 1.5,
     vertex.label.dist = 2.5,
     margin = c(-0.1, -0.1, -0.1, -0.1))


# Graph A ------------------------------------------------------------

# adjacency matrix values
adj_mat_B <- matrix(c(0,1,0,0,0,
                      1,0,1,1,0,
                      0,1,0,1,0,
                      0,1,1,0,0,
                      0,0,0,0,0),
                    nrow = 5,
                    byrow = TRUE) 

# add node names
dimnames(adj_mat_B) = list(c("A", "B", "C", "D", "E"), 
                           c("A", "B", "C", "D", "E")) 

p2 <- graph_from_adjacency_matrix(adj_mat_B,
                                  mode = "undirected",
                                  weighted = TRUE)

plot(p2, layout = layout.fruchterman.reingold,
     vertex.color = "#80b1d3",
     edge.color = "black",
     edge.width = 1.5,
     vertex.label.dist = 2.5,
     margin = c(-0.1, -0.1, -0.1, -0.1))


# Graph C ------------------------------------------------------------

# adjacency matrix values
adj_mat_C <- matrix(c(1,1,0,0,0,
                      1,1,1,1,0,
                      0,1,1,1,0,
                      0,1,1,1,0,
                      0,0,0,0,0),
                    nrow = 5,
                    byrow = TRUE) 

# add node names
dimnames(adj_mat_C) = list(c("A", "B", "C", "D", "E"), 
                           c("A", "B", "C", "D", "E")) 

p3 <- graph_from_adjacency_matrix(adj_mat_C,
                                  mode = "undirected",
                                  weighted = TRUE)


coords <- layout_in_circle(p3)


plot(p3, layout = coords,
     vertex.color = "#80b1d3",
     edge.color = "black",
     edge.width = 1.5,
     vertex.label.dist = 0,
     margin = c(-0.1, -0.1, -0.1, -0.1))

# saving plots
pdf(file = "network_traits_example.pdf", width = 9, height = 3.25, useDingbats=F)
par(mfrow = c(1,3), mar = c(6.1, 2.1, 3.1, 2.1))
plot(p1, layout = coords,
     vertex.color = "#80b1d3",
     edge.color = "black",
     edge.width = 1.5,
     vertex.label.dist = 0,
     margin = c(0, 0, 0, 0.25))
text(-1, 1.5, "A", cex = 1.75, font = 2)
text(0, -1.75, expression(frac(E[L], V) == frac(4,5)), cex = 1.5)
plot(p2, layout = coords,
     vertex.color = "#80b1d3",
     edge.color = "black",
     edge.width = 1.5,
     vertex.label.dist = 0,
     margin = c(0, 0, 0, 0.25))
text(-1, 1.5, "B", cex = 1.75, font = 2)
text(0, -1.75, expression(frac(2*(E[D]), V*(V-1)) == frac(8, 20)), cex = 1.5)
plot(p3, layout = coords,
     vertex.color = "#80b1d3",
     edge.color = "black",
     edge.width = 1.5,
     vertex.label.dist = 0,
     margin = c(0, 0, 0, 0.25))
text(-1, 1.5, "C", cex = 1.75, font = 2)
text(0, -1.75, expression(frac(2*(E[D]) + E[L], V*(V-1) + V) == frac(12, 25)), cex = 1.5)
dev.off()

# calculations ------------------------------------------------------

# (A) host breadth  
ecount(p1)/vcount(p1) # edges are loops only here!

# (B) cross species transmission links (i.e. graph density without loops)
2*(ecount(p2))/(vcount(p2)*(vcount(p2)-1))

# (C) both (i.e. graph density with loops)
(2*(ecount(p2)) + 4) / (vcount(p2)*(vcount(p2)-1) + vcount(p2))

