\name{DRrank}
\alias{DRrank}
\title{ Ranking Regulators by Target Enrichment Density (TED) and Targets' DCL Density (TDD)}
\description{
  The algorithm to rank candidate regulators
}
\usage{
DRrank( DCGs, DCLs, 
	tf2target, expGenes,
	rank.method=c('TED','TDD')[1],
	Nperm=0 )
}
\arguments{
  \item{DCGs}{ a data frame or matrix for DCGs list. }
  \item{DCLs}{ a data frame or matrix for DCLs list. }
  \item{tf2target}{ a data frame or matrix for TF-to-target interaction pairs. }
  \item{expGenes}{ a list for measured genes by array}
  \item{rank.method}{ a character string indicating which ranking method to be utilized. The default is 'TED'. }
  \item{Nperm}{permutation times. If Nperm>0, the permutation step will be implemented for TED and TDD methods. The default value for Nperm is 0.}
}
\details{
DRrank is implemented for ranking potential TFs in terms of their relevance to the phenotypic change or biophysical process of interest. It 
contains two methods: TED, and TDD.

TED, short for 'Target Enrichment Density', employs Binomial Probability model to quantify the enrichment of a TF's targets in the DCG set, 
and as such to evaluate which regulators are more likely to be subject-relevant or even causal. Suppose we sift K DCGs from expression profile 
which contains N genes. If TFi has Ti targets in regulation knowledge, there should be Ti * K / N DCGs appeared in TFi targets 
list randomly. Actually, it is found that TI DCGs are included in TFi's targets list. The larger TI than Ti * K / N is, the more 
targets of TFi enriched, the more likely TFi is a relevant or causative regulator.


TDD, short for 'Targets' DCL Density', uses Clustering Coefficient to quantify the density of DCLs among a regulator's targets, and so to judge 
the importance of a TF. Suppose that TFi has n targets, and that there are k DCLs among these targets. A larger k means more DCLs are 
bridged by the common TFi. We intuitively assume that, if a TF bridged more TF_bridged_DCL it is of more importance (even if the regulator 
is not a DCG).

}
\value{
  A matrix to display TED or TDD scores and ranks.
}
\author{Jing Yang, Hui Yu}

\keyword{ Regulator Rank }

\examples{
data(exprs)
## divide exprs into two parts corresponding to condition 1 
## (exprs.1) and condition 2 (exprs.2) respectively
expGenes<-rownames(exprs)
exprs<-exprs[1:100,]
exprs.1<-exprs[1:100,1:16]
exprs.2<-exprs[1:100,17:63]

data(tf2target)
DCp.res<-DCp(exprs.1,exprs.2,
	link.method = 'qth',cutoff=0.25)
DCe.res<-DCe(exprs.1,exprs.2,
	link.method = 'qth',cutoff=0.25,nbins=10,p=0.1)
DCsum.res<-DCsum(DCp.res,DCe.res,DCpcutoff=0.25,DCecutoff=0.4)

## rank all the potential TFs
data(tf2target)
DRrank.TED.res<-DRrank(DCsum.res$DCGs, DCsum.res$DCLs, 
	tf2target, expGenes,
	rank.method=c('TED','TDD')[1],
	Nperm=0)

DRrank.TED.res[1:3,]

DRrank.TDD.res<-DRrank(DCsum.res$DCGs, DCsum.res$DCLs, 
	tf2target, expGenes,
	rank.method=c('TED','TDD')[2],
	Nperm=0)

DRrank.TDD.res[1:3,]
}
