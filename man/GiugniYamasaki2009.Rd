\name{GiugniYamasaki2009}
\alias{GiugniYamasaki2009}
\docType{data}
\title{
The Policy Impact of Social Movements
}
\description{
This data is used to study the policy impact of social movements.
}
\usage{data(GiugniYamasaki2009)}
\format{
  A data frame with 30 observations on the following 6 variables.
  \describe{
    \item{\code{Case}}{a character vector}
    \item{\code{Protest}}{a numeric vector}
    \item{\code{Opinion}}{a numeric vector}
    \item{\code{Allies}}{a numeric vector}
    \item{\code{Parties}}{a numeric vector}
    \item{\code{Alliances}}{a numeric vector}
    \item{\code{PolicyChange}}{a numeric vector}
  }
}
\details{
This data was collected from a previous study on the policy impact of
antinuclear,ecology, and peace movements in three countries. Giugni and
Sakura replicated the analysis with QCA.

The example section also demonstrates the advantages of using R to do
QCA, because it has strong capacity of data manipulation.
}
\references{
Giugni, Marco and Sakura Yamasaki. 2009. "The Policy Impact of Social Movements: A Replication Through Qualitative Comparative Analysis." Mobilization: An International Quarterly 14 (4) : 467 - 484.
}
\examples{
## ecology movement
ecology <- GiugniYamasaki2009[grep("-E",GiugniYamasaki2009$Case),]
ecology$opinion <- as.numeric(ecology$Opinion>63)
ecology$alliances <- as.numeric(ecology$Alliances>127)
cs_truthTable(ecology,"PolicyChange",c('Protest','opinion','alliances'),case="Case")
reduce(PolicyChange~Protest+opinion+alliances,ecology,case="Case")

## anti-nuclear movement
nuclear <- GiugniYamasaki2009[grep("-N",GiugniYamasaki2009$Case),]
nuclear$opinion <- as.numeric(nuclear$Opinion>51)
nuclear$alliances <- as.numeric(nuclear$Alliances>127)
reduce(PolicyChange~Protest+opinion+alliances,nuclear,case="Case")

## Peace movement
peace <- GiugniYamasaki2009[grep("-P",GiugniYamasaki2009$Case),]
peace$opinion <- as.numeric(peace$Opinion>40)
peace$alliances <- as.numeric(peace$Alliances>127)
reduce(PolicyChange~Protest+opinion+alliances,peace,case="Case")
}
\keyword{datasets}
