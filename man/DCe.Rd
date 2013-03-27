\name{DCe}
\alias{DCe}
\title{ Identify DCGs (Differential Coexpressed Genes) and DCLs (Differential Coexpressed Links) }
\description{
  The algorithm first determines DCLs using a LFC model originally proposed for differential expression analysis, and then determines DCGs 
  with their surrounding links enriched for DCLs.   
}
\usage{
DCe(exprs.1, exprs.2, 
	link.method=c('qth','rth','percent')[1],
	cutoff, 
	r.method, 
	q.method, 
	nbins=20, p=0.1,
	figname=c('LFC.s.jpeg','LFC.d.jpeg'))
}
\arguments{
  \item{exprs.1}{ a data frame or matrix for condition A, with rows as variables (genes) and columns as samples. }
  \item{exprs.2}{ a data frame or matrix for condition B, with rows as variables (genes) and columns as samples. }
  \item{link.method}{ a character string indicating link filtration method, default is 'qth'. }
  \item{cutoff}{ cutoff used for link filtration, can be rth, qth, or percent depending on link.method. must be within [0,1]. }
  \item{r.method}{ a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default) or "spearman", can be abbreviated. }
  \item{q.method}{ a character string indicating which correction method to be utilized. the default is 'BH'.}
  \item{nbins}{ number of x bins for fitting y=a+(b/x). }
  \item{p}{ the cutoff of q-value; must be within [0,1]. }
  \item{figname}{names of figures of the LFC fitting results. }
}
\details{
  DCe is based on the 'Limit Fold Change' or 'LFC' model, a robust statistical method originally proposed for selecting DEGs from microarray 
  data (Mutch et al. 2002). With the analysis units changed from expression values to coexpression values, the LFC method with moderate adaption 
  can be applied to screen for putative DCLs. DCGs with their surrounding links enriched for DCLs are determined through a binomial probability model.
}
\value{
  A list with two components:
  \item{DCGs}{Differential Coexpression Genes}
  \item{DCLs}{Differential Coexpression Links}
}
\references{ 
	Mutch, D.M., Berger, A., Mansourian, R., Rytz, A. and Roberts, M.A. (2002) The limit fold change model: a practical approach for selecting differentially expressed genes from microarray data, BMC Bioinformatics, 3, 17.
  
 Yu H., Liu B-H., et al., Link-specific Quantitative Methods to Identify Differentially Coexpressed Genes and Gene Pairs. Submitted. 2010
 }
\author{ Bao-Hong Liu,  Hui Yu, Jing Yang }

\examples{
data(exprs)
## Identifying DCGs and DCLs by DCe method
DCe.res<-DCe(exprs[1:100,1:16],exprs[1:100,17:63],
	link.method = 'qth',
	cutoff=0.25,
	nbins=20,p=0.1)
DCe.res$DCGs[1:3,]
DCe.res$DCLs[1:3,]
}
\keyword{ DCGs }
\keyword{ DCLs }
