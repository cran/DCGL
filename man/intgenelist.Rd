\name{intgenelist}
\alias{intgenelist}
\docType{data}
\title{ Interesting genes list }
\description{
  Given a list of interesting genes
}
\usage{data(intgenelist)}
\format{
A data frame or matrix with 1 column and several rows of gene symbols.
\itemize{
  \item{intgenelist}{ A data frame with 1 column and several rows. Interesting Gene Symbols. }
}
}
\details{
In the intgenelist, there are 5 Gene symbols in 1 column with the column name 'GeneSymbol' needed.
}
\author{ Jing Yang, Hui Yu }
\examples{
data(intgenelist)
intgenelist # Gene symbol list
colnames(intgenelist) # needed 'GeneSymbol' as the column name
}
\keyword{Interesting Gene Symbol}
