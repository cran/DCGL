\name{expressionBasedfilter}
\alias{expressionBasedfilter}
\title{ Filter genes according to expression level }
\description{
  Genes that have a higher between-sample mean expression signal are retained, while other genes are discarded.
}
\usage{
expressionBasedfilter(exprs)
}
\arguments{
 \item{exprs}{ a data frame or matrix with rows as variables (genes) and columns as samples.}
}
\details{
  Genes which have a Between-Experiment Mean Expression Signal (BEMES) lower than the median of BEMES's of all genes will be filtered out.
}
\value{
A data frame or matrix with a reduced number of rows.
}
\author{ Bao-hong Liu, Hui Yu, Jing Yang }
\examples{
data(exprs)
expressionBasedfilter(exprs)
}
\keyword{ gene filtering }
