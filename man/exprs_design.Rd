\name{exprs_design}
\alias{exprs_design}
\docType{data}
\title{Experiment design of microarray matrix data}
\description{
  A matrix display microarray experiment, which contains rows corresponding to arrays and columns to coefficients to be estimated.
}
\usage{data(exprs_design)}
\format{
A data frame with array samples.
}
\details{
Using numbers to describe microarray experiment design. Usually defined normal samples to be 0, unnormal samples to be 1.
}
\author{ Jing Yang, Hui Yu}
\examples{
data(exprs_design)
exprs_design #define 1 to 16 rows as normal condition, 17 to 63 rows as unnormal condition.
}
\keyword{design}