\name{DCGL-package}
\alias{DCGL-package}
\alias{DCGL}
\docType{package}
\title{
Differential Coexpression Analysis and Differential Regulation Analysis of Gene Expression Microarray Data
}
\description{
DCGL package contains four modules which are Gene filtration module, Link filtration module, Differential CoExpression Analysis (DCEA) 
module and Differential Regulation Analysis (DRA) module. In Gene filtration module, there are expressionBasedfilter and varianceBasedfilter 
functions to filter genes on expression microarray data. In Link filtration module, there are rLinkfilter, percentLinkfilter and qLinkfilter 
functions to filter gene coexpression links in coexpression networks. DCp, DCe, WGCNA, LRC and ASC functions were implemented in DCEA module 
for extracting differentially coexpressed genes (DCGs) and differentially coexpressed links (DCLs). The final step of DCEA module is DCsum to 
determine DCGs and DCLs which come from multiple DCEA methods. In DRA module, there are DRsort, DRplot and DRrank functions to identify 
differentially regulated genes (DRGs) and differentially regulated links (DRLs) and to present some relevant information according to regulation 
knowledge.
}
\details{
\tabular{ll}{
Package: \tab DCGL\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2012-12-14\cr
License: \tab GPL (>2)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Jing Yang, Hui Yu, Bao-Hong Liu

Maintainer: Bao-Hong Liu <bhliu@scbit.org>
}
\references{
Friedrich Leisch,2008 Creating R Packages: A Tutorial
}
\keyword{ package }
\examples{
data(exprs) 
data(tf) ## Transcription Factors (TFs)
data(tf2target) ## TF-to-target relationships
data(exprs_design)

## divide exprs into two parts corresponding to condition 1 (exprs.1) and condition 2 (exprs.2) respectively
exprs.1<-exprs[1:50,1:16] 
exprs.2<-exprs[1:50,17:63]

DCp.res<-DCp(exprs.1,exprs.2,link.method='qth',cutoff=0.25,N=0)
DCe.res<-DCe(exprs.1,exprs.2,link.method='qth',cutoff=0.25,nbins=10,p=0.1)

## combine two Differential CoExpression Analysis results
DCsum.res<-DCsum(DCp.res,DCe.res,DCpcutoff=0.25,DCecutoff=0.4)
DCsum.res$DCGs[1:3,]
DCsum.res$DCLs[1:3,]

## sort out differentially regulated genes and differentially regulated links
DRsort.res<-DRsort(DCsum.res$DCGs,DCsum.res$DCLs,tf2target,exprs[1:50,])

## plot differentially regulated links
DRplot.res<-DRplot(DRsort.res,type='TF_bridged_DCL',vsize=5,asize=0.25,lcex=0.3,ewidth=1,figname=c('TF2target_DCL.pdf','TF_bridged_DCL.pdf'))

## rank regulators by three algorithms
DRrank.res<-DRrank(exprs[1:50,],exprs.1,exprs.2,tf,tf2target,exprs_design,p.value=0.05,DRsort.res)
}
