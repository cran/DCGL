\name{ASC}
\alias{ASC}
\title{
Identify DCGs (Differential Coexpressed Genes) based on 'Average Specific Connection' 
}
\description{
A method to pick out DCGs from microarray data based on 'Average Specific Connection' (ASC) (Choi et al. 2005). 
}
\usage{
ASC(exprs.1, exprs.2, link.method = c("qth", "rth", "percent")[1], cutoff)
}
\arguments{
  \item{exprs.1}{ a data frame or matrix for condition A, with rows as variables (genes) and columns as samples. }
  \item{exprs.2}{ a data frame or matrix for condition B, with rows as variables (genes) and columns as samples. }
  \item{link.method}{ a character string indicating link filtration method, default is 'qth'. }
  \item{cutoff}{ the cutoff of correlation-value, q-value or percent of links after link filtering. must be within [0,1]. }
}
\details{
ASC is the average value of the specific degree of the two conditions. 
}
\value{
	\item{ASC}{ the Average Specific Connections of genes. This measure can be used to rank gene in terms of differential coexpression.}
}
\references{
Choi, J.K., Yu, U., Yoo, O.J. and Kim, S. (2005) Differential coexpression analysis using microarray data and its application to human cancer, Bioinformatics, 21, 4348-4355.
}
\author{
Bao-Hong Liu, Hui Yu
}
\examples{
data(exprs)
ASC(exprs[1:100,1:16],exprs[1:100,17:63],
	link.method = 'qth', cutoff=0.25)
}
\keyword{ DCGs }
