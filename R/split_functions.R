#' create all possible partitions
#'
#' @param n number of elements in set
#' @return list containing partitions
partitions <- function(n){
  if(n > 12){
    stop(sprintf('n = %s requires too many partitions', n))
  }
  part_matrix = partitions::setparts(partitions::restrictedparts(n, 2))
  part_matrix = (part_matrix == 1)
  return(part_matrix[, -1, drop = FALSE])
}

#' return leaf number containing predictions
#'
#' @param rpart.obj rpart fit object
#' @param newdata data.frame with test predictors
#' @return integer vector pointing to leaves containing predictions
predict_nodes <- function(rpart.obj, newdata){
  party_obj = partykit::as.party(rpart.obj)
  pred_nodes = predict(party_obj, newdata = newdata, type = "node")
  return(pred_nodes)
}

#' print ATE at each node
i_treat = function(y, offset, params, wt){
  sfun = function(yval, dev, wt, ylevel, digits){
    paste0(" ATE: ", format(signif(yval, digits)))
  }
  environment(sfun) = .GlobalEnv
  list(y = c(y), parms = NULL, numresp = 1, numy = 1,
       summary = sfun, print = sfun)
}

#' evaluate ATE and number treated
e_treat = function(y, wt, parms) {
  treated = (wt > 0.5)
  wmean = mean(y[treated]) - mean(y[!treated])
  list(label = wmean, deviance = sum(treated))
}

#' split function for max-ate detection
s_treat_max = function(y, wt, x, parms, continuous){
  n = length(y)
  treated = (wt > 0.5)
  if (continuous) {
    goodness = cumulative_te(y, treated, x)
    goodness[is.na(goodness)] = -Inf
    goodness = apply(goodness, 1, max)
    list(goodness = goodness[-n], direction = rep(-1,n-1))
  } else {
    nx = length(unique(x))
    parts = partitions(nx)
    goodness = rep(NA, ncol(parts))
    for(i in 1:ncol(parts)){
      left_g = which(parts[, i])
      right_g = which(!parts[, i])
      left_mean = mean(y[(x %in% left_g) & treated]) -
        mean(y[(x %in% left_g) & !treated])
      right_mean = mean(y[(x %in% right_g) & treated]) -
        mean(y[(x %in% right_g) & !treated])
      goodness[i] = max(left_mean, right_mean)
    }
    max_loc = which.max(goodness)
    max_val = goodness[max_loc]
    best_split_left = which(parts[, max_loc])
    best_split_right = which(!parts[, max_loc])
    direction = c(best_split_left, best_split_right)
    goodness = c(rep(0, length(best_split_left)-1),
                 max_val,
                 rep(0, length(best_split_right)-1))
    names(goodness) = direction[-nx]
    list(goodness=goodness, direction=direction)
  }
}

#' split function for min-ate detection
s_treat_min <- function(y, wt, x, parms, continuous){
  n = length(y)
  treated = (wt > 0.5)
  if (continuous) {
    goodness = cumulative_te(y, treated, x)
    goodness[is.na(goodness)] = Inf
    goodness = apply(goodness, 1, min)
    list(goodness = -goodness[-n], direction = rep(-1,n-1))
  } else {
    nx = length(unique(x))
    parts = partitions(nx)
    goodness = rep(NA, ncol(parts))
    for(i in 1:ncol(parts)){
      left_g = which(parts[, i])
      right_g = which(!parts[, i])

      left_mean = mean(y[(x %in% left_g) & treated]) -
        mean(y[(x %in% left_g) & !treated])

      right_mean = mean(y[(x %in% right_g) & treated]) -
        mean(y[(x %in% right_g) & !treated])

      goodness[i] = min(left_mean, right_mean)
    }

    min_loc = which.min(goodness)
    min_val = goodness[min_loc]
    best_split_left = which(parts[, min_loc])
    best_split_right = which(!parts[, min_loc])
    direction = c(best_split_left, best_split_right)
    goodness = c(rep(0, length(best_split_left)-1),
                 min_val,
                 rep(0, length(best_split_right)-1))
    names(goodness) = direction[-nx]
    list(goodness=-goodness, direction=direction)
  }
}

