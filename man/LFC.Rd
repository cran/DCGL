\name{LFC}
\alias{LFC}
\title{ Select DCLs based on 'Limit Fold Change' model }
\description{
  The limit fold change (LFC) modelis a robust statistical method modeling the relationship between maximum coexpression and log coexpression ratio of genes.
  (Mutch, et al., 2002). The algorithm starts with a set of gene coexpression value pairs each 
  comprising two coexpression values of a gene pair calculated under two different conditions respectively.
}
\usage{
LFC(exprs, nbins, p, sign, figname = "LFC.jpeg")
}
\arguments{
  \item{exprs}{ a two-column data matrix, with column one the coexpression values for one condition and column two those for another. }
  \item{nbins}{ number of x bins for fitting y=a+(b/x). }
  \item{p}{ the fraction at y axis for determining boundary points; must be within [0,1]. }
  \item{sign}{ Specifies the sign type of exprs. Exprs is either 'same-sign', with coexpression pairs of the same sign, or 'different-sign', with coexpression pairs of opposite signs. }
  \item{figname}{ names of figures of the LFC fitting results.}
}
\details{
  According to how the signs of coexpression values are paired, gene links are divided into two parts: the 'same-signed' set and the 'differently-signed' set. From the 'differently-signed' set, the 'correlation-switched' gene links that have two differently-signed coexpression values both surpassing a cutoff value are subtracted, who make the first part of DCLs. The remaining differently-signed gene links in aggregate inherit the title of 'differently-signed' set. 
  For the 'same-signed' set, gene links are binned with respect to their maximum coexpression values, and those links ranked the top p of highest fold changes in each bin are fitted with a simple equation of the form y = a + (b/x); for the 'differently-signed set', the horizontal and vertical axes are exchanged and similar binning and fitting procedures are applied. Links lie above the fitted curves are considered DCLs.
}
\value{
  the identified DCLs will be returned, as a subset of coexpression pairs.
}
\references{ 	
Mutch, D.M., Berger, A., Mansourian, R., Rytz, A. and Roberts, M.A. (2002) The limit fold change model: a practical approach for selecting differentially expressed genes from microarray data, BMC Bioinformatics, 3, 17.
 }
\author{ Bao-Hong Liu, Hui Yu }

\keyword{ DCLs }
