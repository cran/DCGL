\name{exprs}
\alias{exprs}
\docType{data}
\title{Real dataset pulled down from GEO (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17967) }
\description{
  Gene expression dataset, containing 1000 rows and 63 columns.
}
\usage{data(exprs)}
\format{
A data frame with 1000 observations 63 variables.
\itemize{
  \item{exprs}{ A data frame with 1000 observations 63 columns. The expression values. }
}
}
\details{
In the sample gene expression data matrix exprs, it was designed to study gene expression in cirrhotic tissues with (N=16) and without (N=47) HCC.
}
\author{ Jing Yang, Hui Yu }
\examples{
data(exprs)
exprs[,1:16] # exprssion data for condition A
exprs[,17:63] # exprssion data for condition B
row.names(exprs) # gene identifier
}
\keyword{datasets}