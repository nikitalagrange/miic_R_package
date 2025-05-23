% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ScoreNode.R
\name{computeEntropy}
\alias{computeEntropy}
\title{Compute conditional entropy}
\usage{
computeEntropy(X, df_conditioning = NULL)
}
\arguments{
\item{X}{[a vector]
A vector that contains the observational data of the variable.}

\item{df_conditioning}{[a data frame]
The data frame of the observations of the set of conditioning variables}
}
\value{
A list that contains :
\itemize{
\item entro: The estimation of conditional entropy without the
complexity cost.
\item entrok: The estimation of conditional entropy with the
complexity cost (\eqn{entrok = entro + cplx}).}
}
\description{
For discrete variables, the computation is based on the
empirical frequency plus a complexity cost (computed as the
Normalized Maximum Likelihood).
}
\examples{
library(miic)
N <- 1000
X <- sample(x = c(1, 2, 3, 4, 5), size = N, replace = TRUE)
Y <- sample(x = c(1, 2, 3), size = N, replace = TRUE)
res <- computeEntropy(as.factor(X), df_conditioning = as.factor(Y))
message("H(X|Y) = ", res$entro)
message("Hk(X|Y) = ", res$entrok)

}
