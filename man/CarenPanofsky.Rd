\name{CarenPanofsky}
\alias{CarenPanofsky}
\docType{data}
\title{
  Hypothetical data set about TQCA.
}
\description{
  This is a hypothetical data set about TQCA in Caren and Panksky (2005).
}
\usage{data(CarenPanofsky)}
\format{
  A data frame with 17 observations on the following 6 variables.
  \describe{
    \item{\code{public}}{a numeric vector}
    \item{\code{elite}}{a numeric vector}
    \item{\code{affiliate}}{a numeric vector}
    \item{\code{strike}}{a numeric vector}
    \item{\code{e_before_a}}{a numeric vector}
    \item{\code{recognition}}{a numeric vector}
  }
}
\details{
  For condition e_before_a, -9 denotes don't care.
}
\references{
  Caren, Neal and Aaron Panofsky 2005. "TQCA: A Technique for Adding
  Temporality to Qualitative Comparative Analysis." Sociological Methods
  Research 34 (2) : 147-172.
  
  Ragin, Charles C. and Sarah Ilene Strand 2008. "Using Qualitative
  Comparative Analysis to Study Causal Order: Comment on Caren and
  Panofsky (2005)." Sociological Methods Research 36 (4) : 431-441.
}
\examples{
data(CarenPanofsky)
tqca.tt <- cs_truthTable(CarenPanofsky,'recognition',names(CarenPanofsky)[1:5])
reduce(tqca.tt) ## result in Ragin and Strad (2008: 438)
}
\keyword{datasets}
