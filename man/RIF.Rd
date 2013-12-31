\name{RIF}
\alias{RIF}
\title{  Ranking Regulators by Regulator Impact Factor (RIF) Method}
\description{
  The algorithm to rank candidate regulators
}
\usage{
RIF(exprs, exprs.1, exprs.2,
	tf2target, 
	exprs_design,
	p.value)
}
\arguments{
  \item{exprs}{a data frame or matrix for expression dataset, with rows as variables (genes) and columns as samples.}
  \item{exprs.1}{a data frame or matrix for condition A, with rows as variables (genes) and columns as samples.}
  \item{exprs.2}{a data frame or matrix for condition B, with rows as variables (genes) and columns as samples.}
  \item{tf2target}{a data frame or matrix for regulator-to-target interaction pairs.}
  \item{exprs_design}{a data frame or matrix for displaying microarray experiment design.}
  \item{p.value}{ a p value threshold to determine differential expression genes (DEGs).}
}
\details{
RIF method, short for 'Regulator Impact Factor' (Reverter etal. 2010), assesses the change of regulation-accountable expression value of 
Differentially Expressed Genes (DEGs) and correlation coefficient between DEGs and TFs to rank TFs.
}
\value{
  A matrix to display RIF scores and ranks.
}
\references{
Reverter, A., Hudson, N. J., Nagaraj, S. H., Perez-Enciso, M., Dalrymple, B. P. (2010) Regulatory impact factors: unraveling the transcriptional regulation of complex traits from expression data, 26, 896-904.
}
\author{Jing Yang, Hui Yu}

\keyword{ Regulator Rank }

\examples{
data(exprs)
## divide exprs into two parts corresponding to condition 1 
## (exprs.1) and condition 2 (exprs.2) respectively
exprs<-exprs[1:100,]
exprs.1<-exprs[1:100,1:16]
exprs.2<-exprs[1:100,17:63]

DCp.res<-DCp(exprs.1,exprs.2,
	link.method = 'qth',cutoff=0.25)
DCe.res<-DCe(exprs.1,exprs.2,
	link.method = 'qth',cutoff=0.25,nbins=10,p=0.1)
DCsum.res<-DCsum(DCp.res,DCe.res,DCpcutoff=0.25,DCecutoff=0.4)

## rank all the potential TFs
data(tf2target)
data(exprs_design)
\donttest{RIF.res<-RIF(exprs,exprs.1,exprs.2,
	tf2target,
	exprs_design,
	p.value=0.05)
}
\donttest{RIF.res[1:3,]}
}