#' Make a phylogenetic conditioner
#' 
#' Make a phylogenetic conditioner dataset and optionally construct a 
#' formula for one of a selection of R statistical frameworks.
#'
#' @param phy A phylogenetic tree as an object of class [ape::phylo]
#' @param phy_id A one sided formula, where the right hand side should be a vector
#' with ids to match data to the phylogeny. Can be either integer, in which case
#' they are taken to be indexes into `phy$tip.label`, or a character vector
#' in which case it should refer to the tip labels directly. It can refer to a column
#' in the `data` argument (safest) or to any vector in the user's environment,
#' in which case it is assumed to be in the same order as the rows or vectors
#' in `data`, or any other variables referenced in `formula`. You can
#' specify an intercept only model `~1` (the default) which tells `ct_conditioner`
#' to assume any data is already in the same order as the phylogeny's tips.
#' @param formula An optional formula to be augmented with the phylogenetic 
#' conditioning factor
#' @param data An optional data.frame, list or environment containing the variables
#' in `formula` (or `phy_id`)
#' @param framework The framework to construct a `formula` for use with. Ignored
#' if `formula = NULL`. Currently available frameworks are `"lme4"`, `"mgcv"`,
#' `"brms"`, `"inla"`, and the special option `"fixed"`
#'
#' @return A `data.frame` with data that can be used to condition on a 
#' phylogeny in any statistical framework that supports 'mixed' or
#' 'hierarchical' effects.
#' @export
#'
#' @examples
ct_conditioner <- function(phy, phy_id = ~1, 
                           formula = NULL, 
                           data = NULL, 
                           framework = c("lme4", "mgcv", "brms", "inla", "fixed"),
                           helpful = TRUE) {
  
  framework <- match.arg(framework)
  
  phy_id <- model.frame(phy_id, data = data)
  if(length(phy_id) == 0) {
    phy_id <- dplyr::tibble(phy_id = phy$tip.label)
  } 
  
  pairs <- make_pairs(phy)
  
  weights <- make_weights(phy, pairs)
  
  ct_data <- make_data(phy$tip.label, pairs, weights)
  ct_data <- dplyr::left_join(phy_id,
                              ct_data,
                              by = "phy_id")
  ct_data <- as.data.frame(ct_data)
  
  form <- make_formula(framework)
  
  if(is.null(formula) & helpful) {
    message("To use this data with ",
            framework,
            " put this factor on the left hand side of your formula: \n",
            form)
  }
  
  if(is.null(formula)) {
    res <- ct_data
  } else {
    res <- list(data = ct_data, form = form)
  }
  
}

make_pairs <- function(phy) {
  
  depths <- ape::node.depth.edgelength(phy)
  done_tips <- rep(0, ape::Ntip(phy))
  
  sister_pairs <- find_sister_pairs(phy)
  
  ig <- phylo_to_igraph(phy)
  
  node_depths <- ape::node.depth(phy, 2)
  tip_depths <- node_depths[phy$edge[phy$edge[ , 2] %in% seq_len(ape::Ntip(phy)) , 1]]
  
  not_paired <- setdiff(seq_along(phy$tip.label), unlist(sister_pairs))
  
  not_paired <- not_paired[order(tip_depths[not_paired])]
  
  new_pairs <- find_other_pairs(ig, phy, not_paired, not_paired)
  
  list(sisters = sister_pairs, cousins = new_pairs)
  
}

make_weights <- function(phy, pairs) {
  
  cousins <- pairs$cousins
  singletons <- sapply(cousins, length) == 1
  through_root <- sapply(cousins[!singletons], function(x) phytools::fastMRCA(phy, phy$tip.label[x[1]], phy$tip.label[x[2]]))
  through_root <- through_root == ape::Ntip(phy) + 1
  
  zeros <- rep(FALSE, length(cousins))
  zeros[singletons] <- TRUE
  zeros[!singletons][through_root] <- TRUE
  
  cuz_weights <- rep(1, length(cousins))
  cuz_weights[zeros] <- 0
  
  sis_weights <- rep(1, length(pairs$sisters))
  
  list(sisters = sis_weights, cousins = cuz_weights)
}

make_data <- function(tip_labels, pairs, weights) {
  
  ct_sis_dfs <- mapply(make_pair_df, 
                       pairs$sisters, weights$sisters, 
                       MoreArgs = list(tip_labels = tip_labels),
                       SIMPLIFY = FALSE)
  
  ct_sis_dfs <- do.call(rbind, ct_sis_dfs)
  
  ct_sis_dfs$pair_type <- "sisters"
  
  ct_cuz_dfs <- mapply(make_pair_df, 
                       pairs$cousins, weights$cousins, 
                       MoreArgs = list(tip_labels = tip_labels),
                       SIMPLIFY = FALSE)
  
  ct_cuz_dfs <- do.call(rbind, ct_cuz_dfs)
  
  ct_cuz_dfs$pair_type <- "cousins"
  
  rbind(ct_sis_dfs, ct_cuz_dfs)
  
}

make_pair_df <- function(pair, weight, tip_labels) {
  namer <- tip_labels[pair]
  namec <- paste(namer, collapse = "_")
  data.frame(phy_id = namer, ct_group = namec, ct_weights = weight)
}

make_formula <- function(framework = c("lme4", "mgcv", "brms", "inla", "fixed")) {
  
  framework <- match.arg(framework)
  
  form <- switch(framework,
                 lme4 = "(ct_weights|ct_group)",
                 inla = 'f(ct_group, ct_weights, model = "iid")')
  
  form
  
}