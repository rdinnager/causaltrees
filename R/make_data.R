ct_conditioner <- function(phy, form = NULL, data = NULL, framework = c("lme4", "mgcv", "brms", "inla")) {
  tips <- seq_along(phy$tip.label)
  n_tip <- length(tips)
  pairs <- fastmap::faststack()
  while(length(tips) > 0) {
    sisters <- phytools::getSisters(phy, tips[1])
    sisters <- sisters[sisters <= n_tip]
    if(length(sisters) > 0) {
      sisters <- c(tips[1], sisters)
      pairs$push(list(sisters))
      tips <- tips[-sisters]
    }
  }
}