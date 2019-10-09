context('test subgroup tree')
library(data.table)

create_example <- function(n, type = c('binary', 'continuous'),
                           ate = 0.5, seed = 1234){
  type = match.arg(type)
  set.seed(seed)
  df = data.table(group1 = sample(c('A', 'B', 'C', 'D'), n, replace=TRUE),
                  x1 = runif(n),
                  group2 = sample(c('a', 'b', 'c', 'd'), n, replace=TRUE),
                  x2 = runif(n),
                  treated = sample(c(FALSE, TRUE), n, replace=TRUE))
  group_ix_high = with(df, group1 %in% c('A', 'B') & x1 > 0.5)
  group_ix_low = with(df, group1 %in% c('C', 'D') & x1 <= 0.5)
  ix_high_treated = group_ix_high & df$treated
  ix_low_treated = group_ix_low & df$treated

  if(type == 'binary'){
    df$response = rbinom(n, 1, 0.5)
    df$response[ix_high_treated] = rbinom(sum(ix_high_treated), 1, 0.5 + ate/2)
    df$response[ix_low_treated] = rbinom(sum(ix_low_treated), 1, 0.5 - ate/2)
  }else{
    df$response = rnorm(n, 0, 1)
    df$response[ix_high_treated] = rnorm(sum(ix_high_treated), 0.5 + ate/2, 1)
    df$response[ix_low_treated] = rnorm(sum(ix_low_treated), 0.5 - ate/2, 1)
  }
  return(df)
}

test_that('subgroup_tree works',{

  # binary
  df = create_example(n = 1500, type = 'binary', ate = 0.4)
  max_tree = subgroup_tree(df$response, df$treated, df[,1:4],
                           direction = 'max', maxdepth = 2, minbucket=100)
  min_tree = subgroup_tree(df$response, df$treated, df[,1:4],
                           direction = 'min', maxdepth = 2, minbucket=100)
  # continuous
  df = create_example(n = 1500, type = 'continuous', ate = 1)
  max_tree = subgroup_tree(df$response, df$treated, df[,1:4],
                           direction = 'max', maxdepth = 2, minbucket=100)
  min_tree = subgroup_tree(df$response, df$treated, df[,1:4],
                           direction = 'min', maxdepth = 2, minbucket=100)
  df[, mean(response[treated]) - mean(response[!treated]), by = group2]
})





