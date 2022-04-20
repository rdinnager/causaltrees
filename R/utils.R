phylo_to_igraph <- function(phy) {
  temp_phy <- phy
  temp_phy$tip.label <- as.character(seq_along(temp_phy$tip.label))
  temp_phy$node.label <- as.character(length(temp_phy$tip.label) +
                                        seq_len(temp_phy$Nnode))


  ig <- igraph::as.igraph(temp_phy, directed = TRUE)
  ig
}