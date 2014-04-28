\name{CressSnow}
\alias{CressSnow}
\docType{data}
\title{The Outcomes of Homeless Mobilization}
\description{
  The Influence of organization, disruption, political mediation, and
  and framing on the outcome of homeless soical movement organizations.
}
\usage{data(CressSnow)}
\format{
  A data frame with 14 observations on the following 11 variables.
  \describe{
    \item{\code{Caseid}}{The name of Social movement organizations.}
    \item{\code{Viable}}{SMO visibility.}
    \item{\code{Disrupt}}{Disruptive tactics.}
    \item{\code{Allies}}{Sympathetic Allies.}
    \item{\code{City}}{City Support.}
    \item{\code{Diag}}{Diagnostic framming.}
    \item{\code{Prog}}{Prognostic framming.}
    \item{\code{Representation}}{Outcome of representation.}
    \item{\code{Resources}}{Outcome of resources.}
    \item{\code{Rights}}{Outcome of rights.}
    \item{\code{Relief}}{Outcome of relief.}
  }
}

\references{
Cress, Daniel M. and David A. Snow. 2000. "The Outcomes of Homeless
Mobilization: The Influence of Organization, Disruption, Political
Mediation, and and Framing." American Journal of Sociology 105 (4) :
1063-1104.
}
\examples{
cond <- c("Viable", "Disrupt","Allies","City","Diag","Prog")
## see table 5 in p1083
reduce(CressSnow,"Representation",cond,exp="positive",contr="positive",case="Caseid")
reduce(CressSnow,"Resources",cond,exp="positive",contr="positive",case="Caseid")
reduce(CressSnow,"Rights",cond,exp="positive",case="Caseid")
reduce(CressSnow,"Relief",cond,exp="positive",contr="positive",case="Caseid")
## define impact
CressSnow$Impact <- as.numeric(with(CressSnow,Representation+Rights+Relief)>=2)
reduce(CressSnow,"Impact",cond,exp="positive",case="Caseid")
}
\keyword{datasets}
