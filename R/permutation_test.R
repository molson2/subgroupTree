#' permutation test for depth 1 trees
#'
#' Perform a permutation test for a null-effect for the largest ATE discovered
#' from a search using a depth-1 subgroup tree.
#'
#' @param response numeric outcome of interest
#' @param treated boolean vector of treatment assignments
#' @param X data.frame of predictors; must be either numeric or factor types
#' @param direction "max" to search for largest treatment effect, "min" to
#'        search for smallest treatment effect
#' @param n_perm number of permutations to calculate
#' @param ate null average treatment effect to test against
#' @param minbucket minumum number of observations required to consider a split
#' @export
subgroup_perm_test <- function(response, treated, X,
                               direction = c('max', 'min'), n_perm = 1000,
                               ate = 0.0, minbucket = 100){
  direction = match.arg(direction)
  if(!is.numeric(response)){
    stop('response must be numeric')
  }
  agg_op = ifelse(direction == 'max', max, min)

  perms = foreach::`%dopar%`(foreach::foreach(i = seq(n_perm)), {
    treated_perm = sample(treated)
    tscores = rep(NA, ncol(X))
    for(j in 1:ncol(X)){
      tscores[j] = best_split(response, treated_perm, X[, j], agg_op,
                             ate = ate, minbucket = minbucket)
    }
    if(i %% 50 == 0){
      cat('Done with permutation ', i, ' of ', n_perm, '\n')
    }
    agg_op(tscores, na.rm = TRUE)
  })
  unlist(perms)
}

#' compute a t-score for the ATE on a subset of the data defined by "cond"
tscore <- function(response, treated, cond, ate = 0.0){
  nt = sum(cond & treated)
  nc = sum(cond & !treated)
  mt = mean(response[cond & treated])
  mc = mean(response[cond & !treated])
  vart = (mean(response[cond & treated]^2) - mt^2) / nt
  varc = (mean(response[cond & !treated]^2) - mc^2) / nc
  t_score = (mt - mc - ate) / sqrt(vart + varc)
  return(t_score)
}

#' Function to calculate best split
best_split <- function(response, treated, x, agg_op, ate = 0.0,
                       minbucket = 100){
    left_t = NA
    right_t = NA
    agg_t = NA
    if(!is.factor(x)){
      sort_ix = order(x)
      node_scores = cumulative_tscore(response[sort_ix], treated[sort_ix],
                                      x[sort_ix], ate, minbucket)
      agg_t = agg_op(node_scores, na.rm = TRUE)
    }else{
      parts = partitions(nlevels(x))
      for(i in 1:ncol(parts)){
        lgroup = levels(x)[parts[, i]]
        rgroup = levels(x)[!parts[, i]]
        if(sum(x %in% lgroup) >= minbucket){
          left_t = tscore(response, treated, x %in% lgroup, ate)
        }
        if(sum(x %in% rgroup) >= minbucket){
          right_t = tscore(response, treated, x %in% rgroup, ate)
        }
        agg_t = agg_op(agg_t, left_t, right_t, na.rm = TRUE)
      }
    }
    return(agg_t)
}



