\name{DRsort}
\alias{DRsort}
\title{
Identify DRGs (Differential regulated Genes) and DRLs (Differential regulated Links)
}
\description{
  The algorithm is to determine DRGs and DRLs from DCGs and DCLs by TF-to-target interaction knowledge.   
}
\usage{
DRsort(DCGs, DCLs, tf2target, expGenes)
}
\arguments{
  \item{DCGs}{ a data frame or matrix for DCGs list. }
  \item{DCLs}{ a data frame or matrix for DCLs list. }
  \item{tf2target}{ a data frame or matrix for TF-to-target interaction pairs. }
  \item{expGenes}{ a list for measured genes by array }
}
\details{
DRsort, is aimed to sift DCGs and DCLs according to regulation knowledge.

If a DCG is a TF, it is intuitively speculated that its related differential coexpression may be attributed to the change of its regulation 
relationships with its targets. So this type of DCGs are termed Differential Regulation Genes (DRGs). Besides if the upstream TFs of a DCG is 
identified, that DCG is possibly a differentially regulated target of an implicated regulator, and so such DCGs are also kept in the set of DRGs.

If a DCL happens to be a TF-to-target relation, we highlight this DCL because it is the direct attribution to differential regulation. This type 
of DCLs are termed TF2target_DCL. On the other hand, if there are one or more common TFs regulating the two genes of a DCL, we also give priority 
to this DCL because the change in the expression correlation of the two genes could be attributed to the disruption of their co-regulation by the 
common TFs. This type of DCLs are termed TF_bridged_DCLs. TF2target_DCL and TF_bridged_DCL, therefore, together form the set of Differentially 
Regulated Links (DRLs).
}
\value{
  A list with four components:
  \item{DCGs}{Displaying all of Differentially Coexpressed Genes with annotated regulator information whether it is available.}
  \item{DCLs}{Displaying all of Differentially Coexpressed Links with annotated regulator information whether it is available.}
  \item{DRGs}{Differentially Regulated Genes by annotating regulator information. }
  \item{DRLs}{Differentially Regulated Links by annotating regulation pairs information}
  \item{DCG2TF}{Displaying DCGs and upstream TF of DCGs by pairs. It is another format of DRGs. }
  \item{TF_bridged_DCL}{Displaying another format of TF_bridged_DCL for the ease of following investigation.}
}
\author{ Jing Yang, Hui Yu }

\examples{
data(exprs)
data(tf2target)

## Two differential co-expression analysis methods
DCp.res<-DCp(exprs[1:100,1:16],exprs[1:100,17:63],
	link.method = 'qth',cutoff=0.25)
DCe.res<-DCe(exprs[1:100,1:16],exprs[1:100,17:63],
	link.method = 'qth',cutoff=0.25,nbins=20,p=0.1)

## Summarizing DCGs and DCLs from DCp and DCe derived results
DCsum.res<-DCsum(DCp.res,DCe.res,DCpcutoff=0.25,DCecutoff=0.25)

## Sorting out DRGs and DRLs from DCsum-outputted results
\donttest{DRsort.res<-DRsort(DCsum.res$DCGs,DCsum.res$DCLs,tf2target,expGenes)}
## or only sorting out DRGs and DRLs from DCe-outputted results
\donttest{DRsort.res<-DRsort(DCe.res$DCGs,DCe.res$DCLs,tf2target,expGenes)}


## DRGs list
\donttest{DRsort.res$DRGs[1:3,]}

## DRLs list
\donttest{DRsort.res$DRLs[1:3,]}
}
\keyword{ DRGs }
\keyword{ DRLs }
