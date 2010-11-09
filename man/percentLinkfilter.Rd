\name{percentLinkfilter}
\alias{percentLinkfilter}
\title{ To filter gene coexpression links according to the max expression correlation value }
\description{
  Keep a fraction (specified with 'percent') of the links (gene pairs) with the highest max correlation values.
}
\usage{
percentLinkfilter(exprs.1, exprs.2, percent)
}
\arguments{
  \item{exprs.1}{ a data frame or matrix for condition A, with rows as variables (genes) and columns as samples. }
  \item{exprs.2}{ a data frame or matrix for condition B, with rows as variables (genes) and columns as samples. }
  \item{percent}{ fraction of links to be retained.}
}
\details{
  Each gene link is associated with two correlation values (one out of condition A and the other out of condition B) and thus a list of 'maximum absolute values' of the two values is decided. Then these 'maximum absolute values' are sorted in decreasing order. At last, a fraction of gene pairs with the highest max correlation values will be retained.
}
\value{
  A list with two components of data frames, one for filtered data of condition A, the other for the counterpart of condition B.
}
\author{ Bao-Hong Liu, Hui Yu }
\examples{
data(dataC)
percentLinkfilter(dataC[,1:10],dataC[,11:20],percent=0.1)
}
\keyword{ link filtering }
