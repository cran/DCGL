\name{tf2target}
\alias{tf2target}
\docType{data}
\title{The dataset of human Transcription Factors regulate potential target genes }
\description{
  There are 21,4607 TF-to-target interaction pairs downloaded from UCSC (http://genome.ucsc.edu/).
}
\usage{data(tf2target)}
\format{
A data frame with 21,4607 observations 2 variables.
\itemize{
\item{tf2target}{ A data frame with 21,4607 observations 2 columns. TF-to-target interaction pairs. }
}
}
\details{
We downloaded tfbsConsSites and tfbsConsFactors files from UCSC (http://genome.ucsc.edu). tfbsConsSites gives the coordinate information of all TFs 
acted on, tfbsConsFactors gives TFs identifier information. We can predict all the TFs to potential targets relationship based on refGene file which also
came from UCSC and contained genes' position information of hg18. There are 215 TFs and 21,4607 TF-to-target interactions.
}
\author{ Jing Yang, Hui Yu }
\examples{
data(tf2target)
tf2target[1:3,] # TFs (column 1) act on potential target genes (column 2)
}
\keyword{Transcription Factors}
\keyword{Target Genes}