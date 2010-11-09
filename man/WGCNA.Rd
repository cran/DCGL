\name{WGCNA}
\alias{WGCNA}
\title{ To identify DCGs (Differentially-Coexpressed genes) based on the 'Weighted Gene Coexpression Network Analysis' }
\description{
  A method to pick out DCGs from microarray data based on 'Weighted Gene Coexpression Network Analysis' (WGCNA) (Mason, MJ. Et al. 2009; van Nas et al. 2009 ). }
\usage{
WGCNA(exprs.1, exprs.2, power,variant)
}
\arguments{
   \item{exprs.1}{ a data frame or matrix for condition A, with rows as variables (genes) and columns as samples. }
   \item{exprs.2}{ a data frame or matrix for condition B, with rows as variables (genes) and columns as samples. }
   \item{power}{the thresholding parameter, an integer >1.}
   \item{variant}{if the variant is 'WGCNA' the original version is evoked; if it is 'DCp', the length-normalized Euclidean distance is adopted to replace the connectivity difference measure. 
 }

}
\details{
 The 'weighted gene coexpression network analysis' (WGCNA) weights links with correlation coefficients and compares the sums of the correlation coefficients of a gene (Mason, et al., 2009; van Nas, et al., 2009). Correlation coefficients are firstly softly thresholded by a 'power'.
 }
\value{
  \item{WGCNA }{score of 'WGCNA' to identify DCGs}
}
\references{   
Mason, M.J., et al. (2009) Signed weighted gene co-expression network analysis of transcriptional regulation in murine embryonic stem cells, BMC Genomics, 10, 327.

van Nas, A., Guhathakurta, D., Wang, S.S., Yehya, N., Horvath, S., Zhang, B., Ingram-Drake, L., Chaudhuri, G., Schadt, E.E., Drake, T.A., Arnold, A.P. and Lusis, A.J. (2009) Elucidating the role of gonadal hormones in sexually dimorphic gene coexpression networks, Endocrinology, 150, 1235-1249. 
 }
\author{ Bao-Hong Liu, Hui Yu }
\examples{
data(dataC)
WGCNA(dataC[1:100,1:10],dataC[1:100,11:20],power=12,variant='WGCNA')
}
\keyword{ DCGs }
