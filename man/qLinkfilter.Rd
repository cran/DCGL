\name{qLinkfilter}
\alias{qLinkfilter}
\title{ To filter gene coexpression links according to the q-values of expression correlation values }
\description{
  Gene links with q-values of coexpression value pairs in either of two conditions higher than the cutoff (qth) are retained, while the coexpression values of other links are set to zero.
}
\usage{
qLinkfilter(exprs.1, exprs.2, qth)
}
\arguments{
  \item{exprs.1}{ a data frame or matrix for condition A, with rows as variables (genes) and columns as samples. }
  \item{exprs.2}{ a data frame or matrix for condition B, with rows as variables (genes) and columns as samples. }
  \item{qth}{the cutoff of q-value; must be within [0,1].}
}
\details{
  For each of the two conditions, the coexpression values are associated with the corresponding p-values (student T-test of the zero nature of a PCC), and these p-values are sorted and transformed to q-values (or formally, false discovery rates). Gene links with q-values of coexpression values in either of two conditions higher than the cutoff (qth) are reserved.
}
\value{
A list with two components of data frames, one for filtered data of condition A, the other for the counterpart of condition B.
}
\author{ Bao-hong Liu, Hui Yu }
\examples{
data(dataC)
qLinkfilter(dataC[,1:10],dataC[,11:20],qth=0.25)
}
\keyword{link filtering}
