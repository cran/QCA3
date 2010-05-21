\name{Yamasaki2009}
\alias{Yamasaki2009}
\docType{data}
\title{
Movement Impact on Nuclear Energy Policy.
}
\description{
This is a data set about Movement Impact on Nuclear Energy Policy (Yamasaki 2009).
}
\usage{data(Yamasaki2009)}
\format{
  A data frame with 10 observations on the following 8 variables.
  \describe{
    \item{\code{country}}{a factor with levels \code{Belgium} \code{Finland} \code{Germanty} \code{Italy} \code{Netherlands} \code{Spain} \code{Sweden} \code{Switzerland} \code{UK}}
    \item{\code{period}}{a factor with levels \code{1975-1980} \code{1979-1984} \code{1982-1987} \code{1984-1989} \code{1985-1990} \code{1989-1994} \code{1993-1998} \code{1997-2002} \code{1998-2003}}
    \item{\code{soc}}{social movement organisation.}
    \item{\code{pol}}{political allies.}
    \item{\code{ind}}{state-industry relationship.}
    \item{\code{are}}{arena shift.}
    \item{\code{foc}}{focusing event.}
    \item{\code{policyChange}}{Major policy change.}
  }
}
\source{
Manually input by HUANG Ronggui according the truthTable in Yamasaki (2009: 490).
}
\references{
Yamasaki, Sakura 2009. "A Boolean Analysis of Movement Impact on Nuclear Energy Policy." Mobilization: An International Quarterly 14 (4) : 485 - 504.
}
\examples{
data(Yamasaki2009)
reduce(policyChange~soc+pol+ind+are+foc,Yamasaki2009,cases="country")
## or construct truthTable first
cstt <- cs_truthTable(Yamasaki2009, "policyChange",c("soc","pol","ind","are","foc"), cases="country")
tb2 <- reduce(cstt) ## result in table 2
summary(tb2) ## coverage of each configuration
}
\keyword{datasets}
