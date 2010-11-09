\name{DCp}
\alias{DCp}
\title{ To identify DCGs (Differentially-Coexpressed genes) based on the 'Differential Coexpression Profile' }
\description{
   A method to pick out DCGs from microarray data based on a novel concept of 'Differential Coexpression Profiles' (DCp) (Yu et al. 2010).
}
\usage{
DCp(exprs.1, exprs.2, method=Linkfilter.methods,cutoff, N=0)

Linkfilter.methods 
# c("rth", "qth", "percent")
}
\arguments{
  \item{exprs.1}{ a data frame or matrix for condition A, with rows as variables (genes) and columns as samples. }
  \item{exprs.2}{ a data frame or matrix for condition B, with rows as variables (genes) and columns as samples. }
  \item{method}{ link filtering method }
  \item{cutoff}{ cutoff used for link filtration, can be rth, qth, or percent depending on Linkfilter.method. must be within [0,1]. }
  \item{N}{ permutation times. If N>0, the permutation step will be implemented. The default value for N is 0.}
}
\details{
DCp starts with a set of gene coexpression value pairs, where each pair is made up with two coexpression 
values of a gene pair calculated under two different condtions. For a particular gene, a 'condition-specific coexpression profile' 
is defined as the vector of the coexpression values that relate to it in one condition, and the two condition-specific 
coexpression profiles of one gene become the two components of the gene's 'differential coexpression profile'. 
A differential coexpression measure (dC) is calculated from the differential coexpression profile as a length-normalized Euclidean Distance.

Then the samples between the two conditions will be disturbed and the samples will be seperated to two conditions.
Calculate dC of this condition. Repeat the above process for N times . Pool all the dC together to form a null distribution of dC. 
The corresponding statistical significance (p-value) is estimated against null statistics. 

}
\value{
A table with DCGs will be listed including 'dC' value and profile 'length' or 'dC' value, profile 'length','p value' and 'FWER' value.
}
\references{ 
Yu H., Liu B-H., et al., Link-specific Quantitative Methods to Identify Differentially Coexpressed Genes and Gene Pairs. Submitted. 2010
 }
\author{ Bao-Hong Liu, Hui Yu }
\examples{
data(dataC)
DCp(dataC[1:100,1:10],dataC[1:100,11:20],method='qth',cutoff=0.25,N=0)
}
\keyword{ DCGs }
\keyword{ DCLs }
