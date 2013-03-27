\name{DCp}
\alias{DCp}
\title{ Identify DCGs (Differential Coexpressed Genes) based on the 'Differential Coexpression Profile' }
\description{
   A method to pick out DCGs from microarray data based on a novel concept of 'Differential Coexpression Profiles' (DCp) (Yu et al. 2010).
}
\usage{
DCp(exprs.1, exprs.2, 
	r.method = c("pearson", "spearman")[1], 
	link.method = c("qth", "rth", "percent")[1], 
	cutoff = 0.25, 
	N = 0, 
	N.type = c("pooled", "gene_by_gene")[1], 
	q.method = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr")[1])
}
\arguments{
  \item{exprs.1}{ a data frame or matrix for condition A, with rows as variables (genes) and columns as samples. }
  \item{exprs.2}{ a data frame or matrix for condition B, with rows as variables (genes) and columns as samples. }
  \item{r.method}{ a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default) or "spearman", can be abbreviated. }
  \item{link.method}{ a character string indicating link filtering method, default is 'qth'. }
  \item{cutoff}{ cutoff used for link filtration, can be rth, qth, or percent depending on link.method. must be within [0,1]. }
  \item{N}{ permutation times. If N>0, the permutation step will be implemented. The default value for N is 0.}
  \item{N.type}{a character string indicating permutation type, default is 'pooled'.}
  \item{q.method}{ a character string indicating which correction method to be utilized. the default is 'BH'.}
}
\details{
DCp starts with a set of gene coexpression value pairs, where each pair is made up with two coexpression 
values of a gene pair calculated under two different conditions. For a particular gene, a 'condition-specific coexpression profile' 
is defined as the vector of the coexpression values that relate to it in one condition, and the two condition-specific 
coexpression profiles of one gene become the two components of the gene's 'differential coexpression profile'. 
A differential coexpression measure (dC) is calculated from the differential coexpression profile as a length-normalized Euclidean Distance.

Then the samples between the two conditions will be disturbed and the samples will be separated to two conditions.
Calculate dC of this condition. Repeat the above process for N times. Pool all the dC together to form a null distribution of dC. 
The corresponding statistical significance (p-value) is estimated against null statistics (short for pooled). Or calculate p-value of a gene
only in this gene's null distribution of dC (short for gene_by_gene).

}
\value{
A table with DCGs will be listed including 'dC' value and profile 'links' or 'dC' value, profile 'links','p.value' and 'q.vaule' value.
}
\references{
Yu H., Liu B-H., et al., Link-specific Quantitative Methods to Identify Differentially Coexpressed Genes and Gene Pairs. Submitted. 2010
}
\author{ Bao-Hong Liu, Hui Yu, Jing Yang }
\examples{
data(exprs)

## calculate differential coexpressed genes by DCp without permutation
DCp(exprs[1:50,1:16],exprs[1:50,17:63],
	link.method='qth',cutoff=0.25,
	N=0)

## calculate differential coexpressed genes by DCp with 100 times permutation
DCp(exprs[1:50,1:16],exprs[1:50,17:63],
	link.method='qth',cutoff=0.25,
	N=100)
}
\keyword{ DCGs }
\keyword{ DCLs }
