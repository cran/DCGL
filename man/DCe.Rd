\name{DCe}
\alias{DCe}
\title{ To identify DCGs (Differentially-Coexpressed Genes) and DCLs (Differentially-Coexpressed Links) }
\description{
  The algorithm first determines DCLs using a LFC model originially proposed for differential expression analysis, then determines DCGs with their sourrounding links enriched for DCLs.   
}
\usage{
DCe(exprs.1, exprs.2, method=Linkfilter.methods,cutoff, nbins=20, p=0.1,figname=c('LFC.s.jpeg','LFC.d.jpeg'))

Linkfilter.methods 
# c("rth", "qth", "percent")
}
\arguments{
  \item{exprs.1}{ a data frame or matrix for condition A, with rows as variables (genes) and columns as samples. }
  \item{exprs.2}{ a data frame or matrix for condition B, with rows as variables (genes) and columns as samples. }
  \item{method}{ link filtering method }
  \item{cutoff}{ cutoff used for link filtration, can be rth, qth, or percent depending on Linkfilter.method. must be within [0,1]. }
  \item{nbins}{ number of x bins for fitting y=a+(b/x). }
  \item{p}{ the cutoff of q-value; must be within [0,1]. }
  \item{figname}{names of figures of the LFC fitting results}
}
\details{
  DCe is based on the 'Limit Fold Change' or 'LFC' model, a robust statistical method originally proposed for selecting DEGs from microarray data (Mutch et al. 2002). With the analysis units changed from expression values to coexpression values, the LFC method with moderate adaption can be applied to screen for putative DCLs. DCGs with their surrounding links enriched for DCLs are determined through a binomial probability model.
}
\value{
  A list with four components:
  \item{DCGs}{Differentially Coexpressed Genes}
  \item{DCL.same}{Differentially Coexpressed Links of the same sign}
  \item{DCL.diff}{Differentially Coexpressed Links of different signs}
  \item{DCL.switched}{Differentially Coexpressed Links of switched value}
}
\references{ 
	Mutch, D.M., Berger, A., Mansourian, R., Rytz, A. and Roberts, M.A. (2002) The limit fold change model: a practical approach for selecting differentially expressed genes from microarray data, BMC Bioinformatics, 3, 17.
  
 Yu H., Liu B-H., et al., Link-specific Quantitative Methods to Identify Differentially Coexpressed Genes and Gene Pairs. Submitted. 2010
 }
\author{ Bao-Hong Liu, Hui Yu }

\examples{
data(dataC)
DCe(dataC[1:100,1:10],dataC[1:100,11:20],method = 'qth',cutoff=0.25,nbins=20,p=0.1)
}

\keyword{ DCGs }
\keyword{ DCLs }
