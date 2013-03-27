\name{varianceBasedfilter}
\alias{varianceBasedfilter}
\title{To filter genes according to expression variability}
\description{
  Those genes not significantly more variable than the median gene are filtered out.
}
\usage{
varianceBasedfilter(exprs,p)
}
\arguments{
  \item{exprs}{a data frame or matrix with rows for variables (genes) and columns for samples.}
  \item{p}{the probability cut-off of the chi-squared model of the gene-specific variance-like statistics.}
}
\details{
  This is an approximate test of the hypothesis that gene has the same variance as the median variance. A statistical significance criterion 
  based on the variance can be used. If the significance criterion is chosen, then the variance of the log-values for each gene is compared to 
  the median of all the variances. The quantity for each gene compared to a percentile of the chi-square distribution with n-1 degrees of freedom. 
  Those genes not significantly more variable than the median gene are filtered out [BRB-ArrayTools Version 3.7].
}
\value{
A data frame or matrix with a reduced number of rows.
}
\references{
Dr. Richard Simon & Amy Peng Lam, BRB-ArrayTools (v3.7) User's Manual: 'Log expression variation filter'. 
}
\author{Bao-hong Liu, Hui Yu, Jing Yang}
\examples{
data(exprs)
varianceBasedfilter(exprs,0.05)
}
\keyword{gene filtering}