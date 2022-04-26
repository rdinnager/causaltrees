plot_pairs <- function(phy, pair, col = "red") {
  path <- ape::nodepath(phy, pair[1], pair[2])
  edges <- sapply(seq_len(length(path) - 1), function(x) c(path[x], path[x + 1]))
  edges2 <- apply(edges, 2, function(x) apply(phy$edge, 1, function(y) setequal(x, y)))
  edges3 <- apply(edges2, 2, which)
  cols <- rep("black", nrow(phy$edge))
  cols[edges3] <- col
  plot(phy, edge.color = cols)
  invisible(edges3)
}