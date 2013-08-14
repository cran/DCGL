\name{rLinkfilter}
\alias{rLinkfilter}
\title{ Filter gene coexpression links according to correlation coefficient value }
\description{
  Keep a fraction of the links (gene pairs) with the higher correlation values than user given threshold (specified with 'cutoff').
}
\usage{
rLinkfilter(exprs.1, exprs.2, cutoff = 0.8, 
	r.method = c("pearson", "spearman")[1])
}
\arguments{
  \item{exprs.1}{ a data frame or matrix for condition A, with rows as variables (genes) and columns as samples. }
  \item{exprs.2}{ a data frame or matrix for condition B, with rows as variables (genes) and columns as samples. }
  \item{cutoff}{ the threshold of correlation coefficient value. Default is 0.8.}
  \item{r.method}{ a character string to indicate which type of correlation coefficient, 'pearson' (default) or 'spearman', will be computed.}
}
\details{
  Each gene link is associated with two correlation values (one out of condition A and the other out of condition B) and thus a list of 'maximum 
  absolute values' of the two values is decided. Then these links have larger correlation coefficient than threshold are kept.
}
\author{ Bao-Hong Liu, Hui Yu, Jing Yang }
\examples{
data(exprs)
rLinkfilter(exprs[1:100,1:16],exprs[1:100,17:63],cutoff=0.8,r.method='pearson')
}
\keyword{ link filtering }
