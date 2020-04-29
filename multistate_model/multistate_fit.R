# stage-specific recruitment across sites 
# set to identity for no cross-site recruitment
h_block <- block_fn(v_ = c(1, 0,
                           0, 1),
                    c_ = sites, r_ = stages)


# create vec-perumtation matrix
k_fn <- function(i_, j_, stages_, sites_) {
  eij_ = matrix(0, nrow = stages_, ncol = sites_)
  eij_[i_, j_] = 1
  return(kronecker(eij_, t(eij_)))
}

ijs <- expand.grid(1:stages, 1:sites)

k_block = Reduce('+', 
                 lapply(1:nrow(ijs), function(x){
                   k_fn(ijs[x,1], ijs[x,2], stages, sites)
                 }))