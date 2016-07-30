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
    \item{\code{ID}}{a factor}
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
