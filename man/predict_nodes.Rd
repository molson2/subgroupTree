% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/split_functions.R
\name{predict_nodes}
\alias{predict_nodes}
\title{return leaf number containing predictions}
\usage{
predict_nodes(rpart.obj, newdata)
}
\arguments{
\item{rpart.obj}{rpart fit object}

\item{newdata}{data.frame with test predictors}
}
\value{
integer vector pointing to leaves containing predictions
}
\description{
return leaf number containing predictions
}
