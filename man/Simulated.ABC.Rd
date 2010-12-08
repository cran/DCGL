\name{Simulated.ABC}
\alias{Simulated.ABC}
\docType{data}
\title{Simulated datasets based on deliberately-perturbed gene regulation networks (Yu, et al., 2010)}
\description{
  Simulated datasets of 1000-gene scaled(Yu, et al., 2010). 
}
\usage{data(dataC)}
\format{
A data frame with 1000 observations  20 variables.
\itemize{
\item{dataA}{A data frame with 1000 observations 20 columns.  The expression values. }
\item{dataB}{A data frame with 1000 observations 20 columns.  The expression values. }
\item{dataC}{A data frame with 1000 observations 20 columns.  The expression values. }
}
}
\details{
  There are three sets of 1000-gene-scaled simulated datasets used in the paper of Yu et al.(2010) for which the expression data were simulated by software SynTReN using networks originated from E.coli regulatory network. For the three datasets, the underlying networks for condition one were identical but were different for condition two. 
}
\references{
  Yu, H., Liu, B.-H., Ye, Z.-Q., Li, Y.-Y. and Li, Y.-X. (2010) Improved Quantitative Differential Coexpression Analysis to Explore Gene Coexpression Properties in Type II Diabetes,prepare to submit.
}
\author{ Bao-hong Liu, Hui Yu }
\examples{
data(dataC)
dataC[,1:10] # exprssion data for condition 1
dataC[,11:20] # exprssion data for condition 2
row.names(dataC) # gene identifier
}
\keyword{datasets}

