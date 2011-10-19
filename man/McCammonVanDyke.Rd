\name{McCammonVanDyke}
\alias{McCammonVanDyke}
\docType{data}
\title{
Social movement coalition
}
\description{
Data set in McCammon and Van Dyke's meta-analysis of studies on social movement coalition
}
\usage{data(McCammonVanDyke)}
\format{
  A data frame with 24 observations on the following 6 variables.
  \describe{
    \item{\code{ideology}}{ideological alignment}
    \item{\code{threats}}{political threats}
    \item{\code{opportunity}}{political opportunities}
    \item{\code{ties}}{prior social ties between social movement organisations}
    \item{\code{resources}}{plentiful resources}
    \item{\code{coalition}}{outcome of coalition formation}
  }
}
%\details{}
\source{
Manually input from table 13.1 of the reference.
}
\references{
McCammon & Van Dyke. 2010. Applying Qualitative Comparative Analysis to Empirical Studies of Social Movement Coalition Formation. in Van Dyke, N. & McCammon, H. J. (ed.) Strategic alliances: coalition building and social movements. Minneapolis : University of Minnesota Press.
}
\examples{
data(McCammonVanDyke)
reduce(McCammonVanDyke,"coalition",c("ideology","threats","opportunity","ties","resources"),cdontcare='negative')
## result in figure 13.2

reduce(McCammonVanDyke,"coalition",c("ideology","threats","opportunity","ties","resources"),cdontcare='remainder')
## result in figure 13.3
}
\keyword{datasets}
