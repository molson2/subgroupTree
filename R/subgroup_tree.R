#' subgroup detection tree
#'
#' Fits a tree designed to aggressively seek out subsets of the data with
#' large (or small) average treatment effects
#'
#' @param response numeric outcome of interest
#' @param treated boolean vector of treatment assignments
#' @param X data.frame of predictors; must be either numeric or factor types
#' @param direction "max" to search for largest treatment effect, "min" to
#'        search for smallest treatment effect
#' @param ... additional arguments to rpart.control, such as maxdepth, etc.
#' @examples
#' \dontrun{
#'
#' set.seed(123)
#' n = 500
#' p = 5
#' treated = sample(c(TRUE, FALSE), n, replace = TRUE)
#' X = as.data.frame(matrix(rnorm(n * p), n))
#' high_ate = X[, 1] > 0.5
#' response = rbinom(n, 1, ifelse(high_ate & treated, 0.9, 0.5))
#' max_tree = subgroup_tree(response, treated, X, 'max', maxdepth = 2,
#'                         minbucket = 50)
#' print(max_tree)
#' }
#' @export
subgroup_tree <- function(response, treated, X, direction = c('max', 'min'),
                         ...){
  direction = match.arg(direction)
  if(direction == 'max'){
    split_funcs = list(eval = e_treat, split = s_treat_max, init = i_treat)
  }else if(direction == 'min'){
    split_funcs = list(eval = e_treat, split = s_treat_min, init = i_treat)
  }
  control = rpart::rpart.control(minxval = 0, cp = -Inf, ...)
  X$response = response
  X$treated = treated
  fit = rpart::rpart(response ~ . - treated, data=X,
                     weights=ifelse(X$treated, 1, 0.5),
                     method=split_funcs, control=control)
  return(fit)
}

#' @useDynLib subgroupTree, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
