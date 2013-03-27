\name{tf}
\alias{tf}
\docType{data}
\title{Human Transcription Factors}
\description{
  There are 215 human Transcription Factors (TFs) downloaded from UCSC (http://genome.ucsc.edu/).
}
\usage{data(tf)}
\format{
A data frame with 215 observations 1 variables.
\itemize{
\item{tf}{ A list of human TFs. }
}
}
\details{
We downloaded tfbsConsSites and tfbsConsFactors files from UCSC (http://genome.ucsc.edu). tfbsConsSites gives the coordinate information of all TFs 
acted on, tfbsConsFactors gives TFs identifier infomation. We can predict all the TFs to potential targets relationship based on refGene file which also
came from UCSC and contained genes' position information of hg18. There are 215 TFs and 21,4607 TF-to-target interactions.
}
\author{ Jing Yang, Hui Yu }
\examples{
data(tf)
row.names(tf) # TF identifier
}
\keyword{Transcription Factors}

