\name{Hino2009}
\alias{Hino2009}
\docType{data}
\title{
Example data of Time-series QCA.
}
\description{
  Example data in Hino (2009), demonstrating the use of time-serious QCA.
}
\usage{data(Hino2009)}
\format{
  A data frame with 30 observations on the following 5 variables.
  \describe{
    \item{\code{ID}}{a factor with levels \code{A80} \code{A90} \code{B80} \code{B90} \code{D80} \code{D90} \code{FI80} \code{FI90} \code{FR80} \code{FR90} \code{G80} \code{G90} \code{IC80} \code{IC90} \code{IR80} \code{IR90} \code{IT80} \code{IT90} \code{L80} \code{L90} \code{NR80} \code{NR90} \code{NT80} \code{NT90} \code{SWE80} \code{SWE90} \code{SWI80} \code{SWI90} \code{UK80} \code{UK90}}
    \item{\code{VOTE}}{a numeric vector}
    \item{\code{FOREIGN}}{a numeric vector}
    \item{\code{UNEMP}}{a numeric vector}
    \item{\code{DIVERT}}{a numeric vector}
  }
}
\references{
  Hino, Airo 2009. "Time-Series QCA: Studying Temporal Change through Boolean Analysis." Sociological Theory and Methods 24 (2) : 247-265.
}
\examples{
data(Hino2009)
Hino2009$COUNTRY <- gsub("[0-9]*","",Hino2009$ID)
Hino2009$TIME <- gsub("[A-Z]*","",Hino2009$ID)
t2 <- QCA3:::timeDiff(Hino2009,Hino2009$COUNTRY,Hino2009$TIME) ## table 2
t3 <- QCA3:::tsData_tdiff(Hino2009,Hino2009$COUNTRY,Hino2009$TIME)
t3$PRES80s <- 1
t3$PRES80s[match(c("FI","G","SWE","IC","L"),t3$COUNTRY)] <- 0
t3$CONV <- 1 - t3$DIVERT
## manually add PRES80s
rownames(t3) <- t3$ID ## used as case ID.
tt <- cs_truthTable(t3,"VOTE",c("FOREIGN","UNEMP","CONV","PRES80s"))
reduce(tt,expl="positive") ## formula (5)
reduce(tt,expl="positive",remainder="include")## formula (6)
reduce(tt,expl="negative") ## ## formula (7)
reduce(tt,expl="negative",remainder="include") ## ## formula (8)
}
\keyword{datasets}
