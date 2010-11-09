\name{systematicLinkfilter}
\alias{systematicLinkfilter}
\title{A systematic procedure for estimating a cutoff threshold of coexpression networks
 }
\description{
 A procedure that generates a plot of the correlation threshold versus the clustering coefficient for helping define a correlation threshold in coexpression network construction.
}
\usage{
systematicLinkfilter(exprs)
}
\arguments{
 \item{exprs}{ a data frame or matrix with rows as variables (genes) and columns as samples.}
 }
\details{
A systematic procedure for inferring a cutoff threshold of coexpression networks directly from their topological properties. The objective is to automatically select a threshold that preserves as many valid coexpression links as possible, while simultaneously controlling the number of false detections. The procedure is based on comparing the observed clustering coefficient and its randomized counterpart as the number of connections is gradually decreased.
This method is computation-intensive and is optimal for a total gene number less than 500. 
}
\value{
A table of 'correlation threshold' vs. 'clustering coefficient' that may assist the user to determine correlation threshold.
}
\references{
Elo, L.L., Jarvenpaa, H., Oresic, M., Lahesmaa, R. and Aittokallio, T. (2007) Systematic construction of gene coexpression networks with applications to human T helper cell differentiation process, Bioinformatics, 23, 2096-2103.
}
\author{ Bao-hong Liu, Hui Yu }
\examples{
data(dataC)
exprs <- dataC[1:50,]
C_r <- systematicLinkfilter(exprs)
}
\keyword{link filtering}
