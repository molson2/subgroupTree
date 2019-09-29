#include <Rcpp.h>
using namespace Rcpp;

//' calculate cumulative ATE for each potential split
//'
//' @param r response vector
//' @param t indicator variable for treated
//' @param x numeric vector of covariated, sorted low to high
//' @param min_tc smallest number of treated / control units we need
//'        before considering a split
//' @return matrix of left / right daughter ATEs
//' @export
// [[Rcpp::export]]
NumericMatrix cumulative_te(NumericVector r, LogicalVector t,
                            NumericVector x,
                            int min_tc = 1){
  int n = r.size();
  NumericMatrix out(n,2);
  double mean_t = 0.0, mean_c = 0.0;
  double n_t = 0.0, n_c = 0.0;

  // left nodes
  for(int i = 0; i < n; i++){
    if(t[i]){
      mean_t = n_t / (n_t + 1) * mean_t + 1/(n_t + 1) * r[i];
      n_t += 1;
    }else{
      mean_c = n_c / (n_c + 1)* mean_c + 1/(n_c + 1)*r[i];
      n_c += 1;
    }
    if( (n_c < min_tc) | (n_t < min_tc) | (i < n-1 & x[i] == x[i+1]) ){
      out(i,0) = NA_REAL;
    }else{
      out(i,0) = mean_t - mean_c;
    }
  }

  // right nodes
  mean_t = 0.0, mean_c = 0.0;
  n_t = 0.0, n_c = 0.0;
  for(int i = 0; i < (n-1) ; i++){
    if(t[n-i-1]){
      mean_t = n_t / (n_t + 1)* mean_t + 1/(n_t + 1)*r[n-i-1];
      n_t += 1;
    }else{
      mean_c = n_c / (n_c + 1)* mean_c + 1/(n_c + 1)*r[n-i-1];
      n_c += 1;
    }
    if( (n_c < min_tc) | (n_t < min_tc) | (i > 1 & x[n-i-1] == x[n-i-2]) ){
      out(n-i-2,1) = NA_REAL;
    }else{
      out(n-i-2,1) = mean_t - mean_c;
    }
  }
  out(n-1, 1) = NA_REAL;

  return out;
}

//' calculate cumulative ATE t-score for each potential split
//'
//' @param r response vector
//' @param t indicator variable for treatment
//' @param x numeric predictor vector (sorted low to high)
//' @param ate average treatment effect to test
//' @param minbucket only consider splits with at least minbucket obs
//' @return matrix of left / right daughter tscores
//' @export
// [[Rcpp::export]]
NumericMatrix cumulative_tscore(NumericVector r, LogicalVector t,
                                NumericVector x,
                                double ate = 0.0,
                                int minbucket = 100){
  int n = r.size();
  NumericMatrix out(n,2);
  double mean_t = 0.0, mean_c = 0.0;
  double mssq_t = 0.0, mssq_c = 0.0;
  double var_t, var_c, n_t = 0.0, n_c = 0.0;

  // left nodes
  for(int i = 0; i < n; i++){
    if(t[i]){
      mean_t = n_t / (n_t + 1)* mean_t + 1/(n_t + 1)*r[i];
      mssq_t = n_t / (n_t + 1)* mssq_t + 1/(n_t + 1) * r[i] * r[i];
      n_t += 1;
    }else{
      mean_c = n_c / (n_c + 1)* mean_c + 1/(n_c + 1)*r[i];
      mssq_c = n_c / (n_c + 1)* mssq_c + 1/(n_c + 1) * r[i] * r[i];
      n_c += 1;
    }
    if((n_c + n_t < minbucket) | (i < n-1 & x[i] == x[i+1])){
      out(i,0) = NA_REAL;
    }else{
      var_t = (mssq_t - mean_t * mean_t) / n_t;
      var_c = (mssq_c - mean_c * mean_c) / n_c;
      out(i,0) = (mean_t - mean_c - ate) / sqrt(var_t + var_c);
    }
  }

  // right nodes
  mean_t = 0.0, mean_c = 0.0, mssq_t = 0.0, mssq_c = 0.0;
  n_t = 0.0, n_c = 0.0;
  for(int i = 0; i < (n-1) ; i++){
    if(t[n-i-1]){
      mean_t = n_t / (n_t + 1)* mean_t + 1/(n_t + 1)*r[n-i-1];
      mssq_t = n_t / (n_t + 1)* mssq_t + 1/(n_t + 1) * r[n-i-1] * r[n-i-1];
      n_t += 1;
    }else{
      mean_c = n_c / (n_c + 1)* mean_c + 1/(n_c + 1)*r[n-i-1];
      mssq_c = n_c / (n_c + 1)* mssq_c + 1/(n_c + 1) * r[n-i-1] * r[n-i-1];
      n_c += 1;
    }
    if((n_c + n_t < minbucket) | (i > 1 & x[n-i-1] == x[n-i-2])){
      out(n-i-2,1) = NA_REAL;
    }else{
      var_t = (mssq_t - mean_t * mean_t) / n_t;
      var_c = (mssq_c - mean_c * mean_c) / n_c;
      out(n-i-2,1) = (mean_t - mean_c - ate) / sqrt(var_t + var_c);
    }
  }
  out(n-1, 1) = NA_REAL;
  return out;
}

