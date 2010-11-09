\name{DCGL-package}
\alias{DCGL-package}
\alias{DCGL}
\docType{package}
\title{
Differential Coexpression Analysis of Microarray Data 
}
\description{
Functions for basic differential coexpression analyses: gene filtering, link filtering, 
DCG (Differentially-Coexpressed Gene) identification and DCL (Differentially-Coexpressed Links) identification.
Five algorithms,named DCp, DCe, ASC, LRC and WGCNA, are provided for DCG-identification and DCL-identification.
}
\details{
\tabular{ll}{
Package: \tab DCGL\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-08-18\cr
License: \tab GPL (>2)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Bao-Hong Liu, Hui Yu

Maintainer: Bao-hong Liu <bhliu@scbit.org>
}
\references{
Friedrich Leisch,2008 Creating R Packages: A Tutorial
}
\keyword{ package }

\examples{
data(dataC)
exprs.1=dataC[1:100,1:10]
exprs.2=dataC[1:100,11:20]
DCp(exprs.1,exprs.2,method='qth',cutoff=0.25,N=0)
DCe(exprs.1,exprs.2,method='qth',cutoff=0.25,nbins=20,p=0.1)
ASC(exprs.1,exprs.2,method='qth',cutoff=0.25)
LRC(exprs.1,exprs.2,method='qth',cutoff=0.25)
WGCNA(exprs.1,exprs.2,power=12,variant='WGCNA')
}
