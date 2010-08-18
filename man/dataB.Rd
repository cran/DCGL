\name{dataB}
\alias{dataB}
\docType{data}
\title{Simulated dataset based on deliberately-perturbed gene regulation networks (Yu, et al., 2010)}
\description{
  Simulated gene expression dataset, each containing 1000 rows and 20 columns (Yu, et al., 2010).
}
\usage{data(dataB)}
\format{
A data frame with 1000 observations  20 variables.
\itemize{
\item{dataB}{ A data frame with 1000 observations 20 columns.  The expression values. }
}
}
\details{
  We simulated three dataset pairs (denoted dataA, dataB, dataC) using SynTReN based on a predefined E.coli gene regulatory network of a total of 1300 genes (Van den Bulcke, et al., 2006). 
  Specifically, we selected a sub-network of 1000 genes as the original network, 
  and exerted artificial perturbation on 10% of its links as if it was from a different condition. 
  The three groups had different perturbation types. 
  For dataB, we used regulation-elimination (removing a link between a pair of genes). 
  In each data matrix, the first ten columns correspond to one condition and the second ten correspond to the other.
}
\references{
  Yu H., Liu B-H., et al., Link-specific Quantitative Methods to Identify Differentially Coexpressed Genes and Gene Pairs. Submitted. 2010

  Van den Bulcke, T., et al. (2006) SynTReN: a generator of synthetic gene expression data for design and analysis of structure learning algorithms, BMC Bioinformatics, 7, 43.
}
\author{ Bao-hong Liu, Hui Yu }
\examples{
data(dataB)
dataB[,1:10] # exprssion data for condition 1
dataB[,11:20] # exprssion data for condition 2
row.names(dataB) # gene identifier
}
\keyword{datasets}