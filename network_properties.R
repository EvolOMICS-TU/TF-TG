library(igraph)
ls <- list.files(pattern="edgelist.txt")


for(i in ls) {
edge <- read.table(i, header = FALSE, col.names = c("TF_ID", "from", "Gene_ID", "to", "Coord"))
edges <- edge[c(2, 4)]
l1 <- as.data.frame(edges)
graph <- graph_from_data_frame(l1)

# Collect the different centrality measures in a data frame
l2 <- data.frame(degree = degree(graph),
                          closeness = closeness(graph),
                          betweenness = betweenness(graph, directed = TRUE),
                          eigenvalue = eigen_centrality(graph, directed = TRUE)$vector,
                          diameter = diameter(graph, directed = TRUE, unconnected = TRUE, weights = NULL),
                          mean_distance = mean_distance(graph,
                                          weights = NULL,
                                          directed = TRUE,
                                          unconnected = TRUE,
                                          details = FALSE
                                          ),
                          centr_degree = centr_degree(
                                              graph = NULL,
                                              nodes = 0,
                                              mode = "all",
                                              loops = TRUE
                                              ),   
                          transitivity = transitivity(
                          graph,
                          type = "global",
                          vids = NULL,
                          weights = NULL,
                          isolates = "zero"),
                          assortativity_degree = assortativity_degree(
                          graph,
                          directed = TRUE)
                          )
write.table(l2, file = paste(i,"network_properties.csv", sep = "_"), sep = ",")
}                        
