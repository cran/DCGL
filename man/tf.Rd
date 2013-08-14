\name{tf}
\alias{tf}
\docType{data}
\title{Human Transcription Factors}
\description{
  There are 199 human Transcription Factors (TFs) downloaded from UCSC (http://genome.ucsc.edu/).
}
\usage{data(tf)}
\format{
A data frame with 199 observations 1 variables.
\itemize{
\item{tf}{ A list of human TFs. }
}
}
\details{
We downloaded tfbsConsSites and tfbsConsFactors files from UCSC (http://genome.ucsc.edu). tfbsConsSites gives the coordinate information of all TFs 
acted on, tfbsConsFactors gives TFs identifier infomation. We can predict all the TFs to potential targets relationship based on refGene file which also
came from UCSC and contained genes' position information of hg18. There are 199 TFs and 19,9950 TF-to-target interactions.
}
\author{ Jing Yang, Hui Yu }
\examples{
data(tf)
}
\keyword{Transcription Factors}

