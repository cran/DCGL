\name{LRC}
\alias{LRC}
\title{ Identify DCGs (Differential Coexpressed genes) based on 'Log Ratio Connections' }
\description{
  A method to pick out DCGs from microarray data based on 'Log Ratio of Connections' (LRC) (Reverter et al. 2006). 
}
\usage{
LRC(exprs.1, exprs.2, link.method = c("qth", "rth", "percent")[1], cutoff)
}
\arguments{
  \item{exprs.1}{ a data frame or matrix for condition A, with rows as variables (genes) and columns as samples. }
  \item{exprs.2}{ a data frame or matrix for condition B, with rows as variables (genes) and columns as samples. }
  \item{link.method}{ a character string indicating link filtering method, default is 'qth'. }
  \item{cutoff}{ cutoff used for link filtration, can be rth, qth, or percent depending on Linkfilter.method. must be within [0,1]. }
}
\details{
 'Log Ratio of Connections' (LRC) calculates the logarithm of the ratio of the connectivities of a gene between two conditions (Reverter, et al., 2006). A connectivity of zero is changed to one.
}
\value{
  \item{LRC}{the log Ratio Connections of genes. This measure can be used to rank gene in terms of differential coexpression.}
}
\references{   
Reverter, A., Ingham, A., Lehnert, S.A., Tan, S.H., Wang, Y., Ratnakumar, A. and Dalrymple, B.P. (2006) Simultaneous identification of differential gene expression and connectivity in inflammation, adipogenesis and cancer, Bioinformatics, 22, 2396-2404.
}
\author{ Bao-Hong Liu, Hui Yu, Jing Yang }
\examples{
data(exprs)
LRC(exprs[1:100,1:16],exprs[1:100,17:63],link.method = 'qth', cutoff=0.25)
}
\keyword{ DCGs}
