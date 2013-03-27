\name{qLinkfilter}
\alias{qLinkfilter}
\title{ Filter gene coexpression links according to the q-values of expression correlation values }
\description{
  Gene links with q-values of coexpression value pairs in either of two conditions higher than the cutoff are retained, while the coexpression values of other links are set to zero.
}
\usage{
qLinkfilter(exprs.1, exprs.2, cutoff = 0.25, r.method = c("pearson", "spearman")[1], q.method = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr")[1])
}
\arguments{
  \item{exprs.1}{ a data frame or matrix for condition A, with rows as variables (genes) and columns as samples. }
  \item{exprs.2}{ a data frame or matrix for condition B, with rows as variables (genes) and columns as samples. }
  \item{cutoff}{ the cutoff of q-value; must be within [0,1].}
  \item{r.method}{ a character string indicating 'pearson' (default) or 'spearman' correlation coefficient will be computed.}
  \item{q.method}{ a character string indicating adjust method of p-value, default is 'BH'.}
}
\details{
  For each of the two conditions, the coexpression values are associated with the corresponding p-values (student T-test of the zero nature of a PCC), 
  and these p-values are sorted and transformed to q-values (false discovery rates). Gene links with q-values of coexpression values in 
  either of two conditions lower than the cutoff are reserved.
}
\value{
A list with two components of data frames, one for filtered data of condition A, the other for the counterpart of condition B.
}
\author{ Bao-hong Liu, Hui Yu, Jing Yang }
\examples{
data(exprs)
qLinkfilter(exprs[1:100,1:16],exprs[1:100,17:63],cutoff=0.25,r.method='pearson',q.method='BH')
}
\keyword{link filtering}
