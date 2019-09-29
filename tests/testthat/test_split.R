context('test cumulative split functions')

test_that("cumulative ates and tscores work", {

  # problem setup
  set.seed(123)
  n = 100
  x = rnorm(n)
  t = sample(c(FALSE, TRUE), n, replace=T)
  y = rnorm(n, 1, 0.5)

  sort_ix = order(x)

  tscore <- function(y, t, cond){
    mt = mean(y[cond & t])
    mc = mean(y[cond & !t])
    vart = (mean(y[cond & t]^2) - mt^2) / sum(cond & t)
    varc = (mean(y[cond & !t]^2) - mc^2) / sum(cond & !t)
    return((mt - mc) / sqrt(vart + varc))
  }

  ate <- function(y, t, cond){
    mt = mean(y[cond & t])
    mc = mean(y[cond & !t])
    return(mt - mc)
  }

  xvals = sort(unique(x))
  tscores = rep(NA, length(xvals), 1)
  ates = rep(NA, length(xvals), 1)

  for(i in 1:length(xvals)){
    left_t = tscore(y, t, x <=xvals[i])
    right_t = tscore(y, t, x > xvals[i])
    tscores[i] = max(left_t, right_t, na.rm=TRUE)
    left_t = ate(y, t, x <= xvals[i])
    right_t = ate(y, t, x > xvals[i])
    ates[i] = max(left_t, right_t, na.rm = TRUE)
  }

  node_tscores = cumulative_tscore(y[sort_ix], t[sort_ix], x[sort_ix], 0.0, 1)
  node_ates = cumulative_te(y[sort_ix], t[sort_ix], x[sort_ix], 1)

  node_tscores = pmax(node_tscores[,1], node_tscores[,2], na.rm = TRUE)
  node_ates = pmax(node_ates[,1], node_ates[,2], na.rm = TRUE)

  expect_equal(node_tscores, tscores)
  expect_equal(node_ates, ates)

})

