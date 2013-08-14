\name{percentLinkfilter}
\alias{percentLinkfilter}
\title{ Filter gene coexpression links according to the top percent of expression correlation coefficient value }
\description{
  Keep a fraction (specified with 'cutoff') of the links (gene pairs) with the highest max correlation values.
}
\usage{
percentLinkfilter(exprs.1, exprs.2, cutoff = 0.25, 
	r.method = c("pearson", "spearman")[1])
}
\arguments{
  \item{exprs.1}{ a data frame or matrix for condition A, with rows as variables (genes) and columns as samples. }
  \item{exprs.2}{ a data frame or matrix for condition B, with rows as variables (genes) and columns as samples. }
  \item{cutoff}{ fraction of links to be retained. Default is 0.25.}
  \item{r.method}{ a character string to indicating 'pearson' (default) or 'spearman' correlation coefficient will be computed.}
}
\details{
  Each gene link is associated with two correlation values (one out of condition A and the other out of condition B) and thus a list of 'maximum 
  absolute values' of the two values is decided. Then these 'maximum absolute values' are sorted in decreasing order. At last, a fraction of gene 
  pairs with the highest max correlation values will be retained.
}
\value{
  A list with two components of data frames, one for filtered data of condition A, the other for the counterpart of condition B.
}
\author{ Bao-Hong Liu, Hui Yu, Jing Yang }
\examples{
data(exprs)
percentLinkfilter(exprs[1:100,1:16],exprs[1:100,17:63],cutoff=0.25,r.method='pearson')
}
\keyword{ link filtering }
