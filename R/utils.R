phylo_to_igraph <- function(phy) {
  temp_phy <- phy
  temp_phy$tip.label <- as.character(seq_along(temp_phy$tip.label))
  temp_phy$node.label <- as.character(length(temp_phy$tip.label) +
                                        seq_len(temp_phy$Nnode))


  ig <- igraph::as.igraph(temp_phy, directed = TRUE)
  ig
}

find_sister_pairs <- function(phy) {
  tips <- seq_along(phy$tip.label)
  n_tip <- length(tips)
  pairs <- fastmap::faststack()
  done <- rep(0, length(tips))
  for(i in seq_along(tips)) {
    if(i %in% which(done == 0)) {
      sisters <- phytools::getSisters(phy, tips[i])
      sisters <- sisters[sisters <= n_tip]
      if(length(sisters) > 0) {
        sisters <- c(tips[i], sisters)
        pairs$push(list(sisters))
        done[sisters] <- 1
      }
    }
  }
  unlist(pairs$as_list(), recursive = FALSE)
}

find_other_pairs <- function(ig, phy, to_pair, pool) {
  pool_done <- rep(0, length(pool))
  to_pair_done <- rep(0, length(to_pair))
  new_pairs <- fastmap::faststack()
  
  for(i in seq_along(to_pair)) {
    if(i %in% which(to_pair_done == 0)) {
      not_paired2 <- setdiff(pool, pool[which(pool_done == 1)])
      other_tips <- as.character(setdiff(not_paired2, to_pair[i]))
      paths <- igraph::shortest_paths(ig,  as.character(to_pair[i]), 
                                      other_tips,
                                      mode = "all",
                                      output = "vpath")
      path_lens <- sapply(paths$vpath, length)
      if(length(path_lens) > 0) {
        pairing <- as.numeric(other_tips[path_lens == min(path_lens)])
        if(length(pairing) > 1) {
          dists <- sapply(pairing, function(x) phytools::fastDist(phy, phy$tip.label[x], phy$tip.label[to_pair[i]]))
          if(var(dists) != 0) {
            pairing <- pairing[1]
          }
        }
        pair <- c(to_pair[i], pairing)
      } else {
        pair <- to_pair[i]
      }
      
      
      new_pairs$push(list(pair))
      to_pair_done[to_pair %in% pair] <- 1
      pool_done[pool %in% pair] <- 1
    }
  }
  unlist(new_pairs$as_list(), recursive = FALSE)
}